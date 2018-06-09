#include <unistd.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#define ONLINE_BUFFER_LENGTH 7500000
#define ONLINE_SHM_KEY 19999
struct OMBuffer_t {
  unsigned long length;
  unsigned long dead;
  unsigned long pWrite;
  unsigned long pRead;
  unsigned long data[ONLINE_BUFFER_LENGTH];
};

struct buffer_data {
  char *base;
  int size;
  int offset;
};

struct timing_data {
  time_t last_t;
  int count;
  float max_rate;
  int precalib;
};

volatile int end_now = 0;

void sig_int (int s) {
  end_now = 1;
  fprintf (stderr, "signal received, exiting\n");
}

void *attach_shm () {
  int shm_handle;
  void *ptr;
  while (1) {
    shm_handle = shmget (ONLINE_SHM_KEY, 0, 0);
    if (shm_handle < 0) {
      if (errno == ENOENT) usleep (300000);  /* if not yet created sleep ad try again */
      else {
	perror ("shmget:");
        exit (1);
      }
    } 
    else break;
  }

  ptr = shmat (shm_handle, 0, 0);
  if (ptr == (void *)-1) {
    perror ("shmat");
    exit (1);
  }

  return ptr;
}

void copy_event (struct buffer_data *buffer, struct OMBuffer_t *shm, int size) {
  int ev_words = size / 4;
  if (buffer->size - buffer->offset <= size) {
    buffer->size = buffer->size + size * 2;
    buffer->base = realloc (buffer->base, buffer->size);
    if (!buffer->base) {
      fprintf (stderr, "realloc failed\n");
      exit (0);
    }
  }

  if (shm->pRead + ev_words < ONLINE_BUFFER_LENGTH) memcpy (buffer->base + buffer->offset, shm->data + shm->pRead, size);
  else {
    int up_to_end = (ONLINE_BUFFER_LENGTH - shm->pRead) * 4;
    memcpy (buffer->base + buffer->offset, shm->data + shm->pRead, up_to_end);
    memcpy (buffer->base + buffer->offset + up_to_end, shm->data, size - up_to_end);
  }

  buffer->offset += size;
  shm->pRead = (shm->pRead + ev_words) % ONLINE_BUFFER_LENGTH;
}

void flush_and_write (struct buffer_data *buffer, struct OMBuffer_t *shm) {
  while (((shm->pWrite + ONLINE_BUFFER_LENGTH - shm->pRead) % ONLINE_BUFFER_LENGTH) > 5000) {
    int size = shm->data[shm->pRead];
    shm->pRead = (shm->pRead + size / 4) % ONLINE_BUFFER_LENGTH;
  }

  if (buffer->offset) fwrite (buffer->base, buffer->offset, 1, stdout);
  buffer->offset = 0;
}

void write_event (struct OMBuffer_t *shm, int size) {
  int ev_words = size / 4;
  if (shm->pRead + ev_words < ONLINE_BUFFER_LENGTH) fwrite (shm->data + shm->pRead, size, 1, stdout);
  else {
    int up_to_end = (ONLINE_BUFFER_LENGTH - shm->pRead) * 4;
    fwrite (shm->data + shm->pRead, up_to_end, 1, stdout);
    fwrite (shm->data, size - up_to_end, 1, stdout);
  }
  shm->pRead = (shm->pRead + ev_words) % ONLINE_BUFFER_LENGTH;
}

int main (int argc, char *argv[]) {
  struct OMBuffer_t *shm;
  struct buffer_data buffer;
  struct timing_data timing;

  if (argc == 2) {
    timing.max_rate = atof (argv[1]);
    if (timing.max_rate <= 0) {
      fprintf (stderr, "rate must be >= 0 (%s)\n", argv[1]);
      exit (0);
    }
    buffer.size = (int)(50000 * timing.max_rate * 10);
    buffer.base= (char *)malloc (buffer.size);
    if (!buffer.base) {
      fprintf (stderr, "malloc failed\n");
      exit (0);
    }
    buffer.offset = 0;
    timing.precalib = 0;
    timing.last_t = 0;
  } else timing.max_rate = 0;

  shm = attach_shm ();

  if (signal (SIGINT, sig_int) == SIG_ERR) {
    perror ("signal");
    exit (1);
  }


  while (1) {
    int free_words;
    unsigned long ev_size;
#if DEBUG
    unsigned long run_number, ev_number;
#endif

    if (shm->dead) break;
    if (end_now) break;

      /* check not to be too close to write limit */
    free_words = (shm->pWrite + ONLINE_BUFFER_LENGTH - shm->pRead) % ONLINE_BUFFER_LENGTH;
    if (free_words < 5000) {
      usleep (100000);	/* wait for 0.1s */
      continue;
    }

      /* read event by event */
    ev_size = shm->data[shm->pRead];
#if DEBUG
    run_number = shm->data[(shm->pRead + 1) % ONLINE_BUFFER_LENGTH];
    ev_number = shm->data[(shm->pRead + 2) % ONLINE_BUFFER_LENGTH];
    fprintf (stderr, "event %lu size %lu run %lu\n", ev_number, ev_size, run_number);
#endif

    if (timing.max_rate && timing.precalib > 1050) {
      time_t t = time (0);
      if (!timing.last_t) timing.last_t = time (0);
      if (t - timing.last_t < 10) {
	timing.count ++;
	if (timing.count < timing.max_rate * 10) copy_event (&buffer, shm, ev_size);
	else if (timing.count == (int)(timing.max_rate * 10)) flush_and_write (&buffer, shm);
	else shm->pRead = (shm->pRead + ev_size / 4) % ONLINE_BUFFER_LENGTH; // do nothing, just skip event
      } else {
	timing.last_t = t;
	timing.count = 0;
      }
    } else {
      write_event (shm, ev_size);
      timing.precalib ++;
    }
  }

  fclose (stdout);

  sleep (3);
  return 0;
}
