//#include <unistd.h>
//#include <sys/types.h>
//#include <sys/stat.h>
//#include <sys/socket.h>
#include <stdlib.h>
#include <stdio.h>
//#include <assert.h>

#include <QtGui>
#include <QTextBrowser>
//#include <QSocketNotifier>

#include "mainWindow.h"
#include "dbdialog.h"
#include "bchmondialog.h"
#include "lmsmondialog.h"
#include "mchmondialog.h"
#include "mmsmondialog.h"
#include "dmpmondialog.h"
#include "dmpmumondialog.h"
#include "hvmondialog.h"
// #include "chsmondialog.h"
// #include "cbsmondialog.h"
// #include "cpsmondialog.h"
// #include "tmsmondialog.h"

mainWindow *mainWindow::me = 0;
int mainWindow::cur_run    = -1;
FILE *term_stream;

mainWindow::mainWindow()
{
  me = this;

  if (!(term_stream = fopen("session_screen.log", "w")))
    QMessageBox::warning(this, tr("WARNING"), tr("<h2>dbMon::mainWindow:</h2>" "<h2>Can not open session_screen.log   </h2>"));

  // db manager address and the current run's number initialization;
  db_mng = pg_dbi::get();
  if (cur_run < 0){
    int lines, columns;
    char the_query[4096];

    sprintf(the_query, "SELECT \"RunNumber\" FROM \"LabenPrecalibDecodingQuality\" ORDER BY \"RunNumber\" ");
    if (pg_dbi::get()->open_connection(3)){
      pg_dbi::get()->query(3, the_query);
      pg_dbi::get()->query_dimensions(lines, columns);
      cur_run = QString::fromAscii(pg_dbi::get()->query_value(lines-1, 0)).toInt();
      pg_dbi::get()->close_connection(3);
    }

    if (pg_dbi::get()->open_connection(2)){
      while (cur_run > 5003){
	sprintf(the_query, "SELECT \"RunNumber\" FROM \"LabenChannelsProperties\" WHERE \"RunNumber\"=%d ", cur_run);
	pg_dbi::get()->query(2, the_query);
        pg_dbi::get()->query_dimensions(lines, columns);
	if (lines > 0) break;
	else --cur_run;
      }
      pg_dbi::get()->close_connection(2);
    }
  }

  textBrowser = new QTextBrowser;
  //textBrowser->setTextColor(Qt::white); textBrowser->setTextBackgroundColor(Qt::black);
  textBrowser->setStyleSheet("color: white; background-color: black");
  textBrowser->setAutoFormatting(QTextEdit::AutoNone);
  textBrowser->setAlignment(Qt::AlignLeft);
  //textBrowser->setLineWrapMode(QTextEdit::NoWrap);

  setCentralWidget(textBrowser);

  createActions();
  createMenus();
  createContextMenu();
  createToolBars();
  createStatusBar();

  dbMngDialog  = 0;
  bchMonDialog = 0;
  lmsMonDialog = 0;
  mchMonDialog = 0;
  mmsMonDialog = 0;
  dmpMonDialog = 0;
  dmpmuMonDialog = 0;
  hvMonDialog  = 0;
//   chsMonDialog = 0;
//   cbsMonDialog = 0;
//   cpsMonDialog = 0;
//   tmsMonDialog = 0;

  setWindowIcon(QIcon(":/images/icon.xpm"));
  setWindowTitle(tr("BOREXINO DataBase Monitoring tool"));
  resize(640, 480); //adjustSize();
}

mainWindow::~mainWindow()
{
  if (db_mng) delete db_mng;
  if (term_stream) fclose(term_stream);
  me = 0;
}

void mainWindow::dispatchMsg(const QString& str, int qualifier){
  textBrowser->insertPlainText(str);
  //textBrowser->update();
  if (term_stream){
    fprintf(term_stream, "%s", str.toLocal8Bit().constData());
    fflush(term_stream);
  }
}

void mainWindow::createActions()
{

  openAllDBAction = new QAction(tr("Open &All"), this);
  //openAllDBAction->setIcon(QIcon(":/images/icon.png"));
  openAllDBAction->setShortcut(tr("Ctrl+A"));
  openAllDBAction->setStatusTip(tr("Open all databases"));
  openAllDBAction->setData(-1);
  connect(openAllDBAction, SIGNAL(triggered()), this, SLOT(openDB()));

  closeAllDBAction = new QAction(tr("&Quit"), this);
  //closeAllDBAction->setIcon(QIcon(":/images/icon.png"));
  closeAllDBAction->setShortcut(tr("Ctrl+Q"));
  closeAllDBAction->setStatusTip(tr("Close all databases and Quit"));
  closeAllDBAction->setData(-1);
  connect(closeAllDBAction, SIGNAL(triggered()), this, SLOT(closeDB()));

  for (int i = 0; i < NDB; ++i) {
    QString dbName(db_names[i]);

    openDBActions[i]    = new QAction("&Open", this);
    operateDBActions[i] = new QAction("&Manage", this);
    closeDBActions[i]   = new QAction("&Close", this);

    openDBActions[i]->setEnabled(true);
    openDBActions[i]->setStatusTip(tr("Open %1 connection").arg(dbName));
    openDBActions[i]->setData(i);

    operateDBActions[i]->setEnabled(false);
    operateDBActions[i]->setStatusTip(tr("%1 tables' content").arg(dbName));
    operateDBActions[i]->setData(i);

    closeDBActions[i]->setEnabled(false);
    closeDBActions[i]->setStatusTip(tr("Close %1 connection").arg(dbName));
    closeDBActions[i]->setData(i);

    connect(openDBActions[i],    SIGNAL(triggered()), this, SLOT(openDB()));
    connect(operateDBActions[i], SIGNAL(triggered()), this, SLOT(operateDB()));
    connect(closeDBActions[i],   SIGNAL(triggered()), this, SLOT(closeDB()));
  }

  bchMonAction = new QAction(tr("&1 Bad/Off in precalibration"), this);
  connect(bchMonAction, SIGNAL(triggered()), this, SLOT(bchMon()));

  lmsMonAction = new QAction(tr("&2 Laben Channels Status"), this);
  connect(lmsMonAction, SIGNAL(triggered()), this, SLOT(lmsMon()));

  mchMonAction = new QAction(tr("&3 Muon Trigger rates"), this);
  connect(mchMonAction, SIGNAL(triggered()), this, SLOT(mchMon()));

  mmsMonAction = new QAction(tr("&4 Muon Channels Status"), this);
  connect(mmsMonAction, SIGNAL(triggered()), this, SLOT(mmsMon()));

  hvMonAction  = new QAction(tr("&5 HV Channels Status"), this);
  connect(hvMonAction, SIGNAL(triggered()), this, SLOT(hvMon()));

  dmpMonAction = new QAction(tr("&Laben Broken Channels"), this);
  connect(dmpMonAction, SIGNAL(triggered()), this, SLOT(dmpMon()));

  dmpmuMonAction = new QAction(tr("&Muon Broken Channels"), this);
  connect(dmpmuMonAction, SIGNAL(triggered()), this, SLOT(dmmMon()));

//   chsMonAction = new QAction(tr("&5 Laben Charge Status"), this);
//   connect(chsMonAction, SIGNAL(triggered()), this, SLOT(chsMon()));
//
//   cbsMonAction = new QAction(tr("&6 Laben Charge Base Status"), this);
//   connect(cbsMonAction, SIGNAL(triggered()), this, SLOT(cbsMon()));
//
//   cpsMonAction = new QAction(tr("&7 Laben Charge Peak Status"), this);
//   connect(cpsMonAction, SIGNAL(triggered()), this, SLOT(cpsMon()));
//
//   tmsMonAction = new QAction(tr("&8 Laben Timing Status"), this);
//   connect(tmsMonAction, SIGNAL(triggered()), this, SLOT(tmsMon()));

  helpAction = new QAction(tr("&Help"), this);
  helpAction->setStatusTip(tr("Help"));
  connect(helpAction, SIGNAL(triggered()), this, SLOT(help()));

  aboutAction = new QAction(tr("&About"), this);
  aboutAction->setStatusTip(tr("dbMon version"));
  connect(aboutAction, SIGNAL(triggered()), this, SLOT(about()));

  aboutQtAction = new QAction(tr("&Qt specs"), this);
  aboutQtAction->setStatusTip(tr("QT library specifications"));
  connect(aboutQtAction, SIGNAL(triggered()), qApp, SLOT(aboutQt()));
}

void mainWindow::createMenus()
{
  dbMenu = menuBar()->addMenu(tr("Data&Bases"));
  dbMenu->addAction(openAllDBAction);
  dbMenu->addSeparator();
  for (int i = 0; i < NDB; ++i) {
    QString dbName(db_names[i]);
    dbButtons[i] = dbMenu->addMenu(tr("&%1 %2").arg(i+1).arg(dbName));
    dbButtons[i]->addAction(openDBActions[i]);
    dbButtons[i]->addAction(operateDBActions[i]);
    dbButtons[i]->addAction(closeDBActions[i]);
  }
  dbMenu->addSeparator();
  dbMenu->addAction(closeAllDBAction);

  //queryMenu = menuBar()->addMenu(tr("&Query"));

  monitorMenu = menuBar()->addMenu(tr("&Tools"));
  monitorMenu->addAction(bchMonAction);
  monitorMenu->addAction(lmsMonAction);
  monitorMenu->addAction(mchMonAction);
  monitorMenu->addAction(mmsMonAction);
  monitorMenu->addAction(hvMonAction);
  //monitorMenu->addAction(chsMonAction);
  //monitorMenu->addAction(cbsMonAction);
  //monitorMenu->addAction(cpsMonAction);
  //monitorMenu->addAction(tmsMonAction);

  dumpMenu = menuBar()->addMenu(tr("&Dump"));
  dumpMenu->addAction(dmpMonAction);
  dumpMenu->addAction(dmpmuMonAction);

  menuBar()->addSeparator();

  helpMenu = menuBar()->addMenu(tr("&Help"));
  helpMenu->addAction(helpAction);
  helpMenu->addAction(aboutAction);
  helpMenu->addAction(aboutQtAction);
}


void mainWindow::createContextMenu()
{
  //testBrowser->addAction(Action);
  //textBrowser->addAction(Action);
  //textBrowser->addAction(Action);
  textBrowser->setContextMenuPolicy(Qt::ActionsContextMenu);
}

void mainWindow::createToolBars()
{
//   browserToolBar = addToolBar(tr("&Terminal Browser"));
//   browserToolBar->addAction(Action);
//   browserToolBar->addAction(Action);
//   browserToolBar->addAction(Action);
}

void mainWindow::createStatusBar()
{
  statusLabel = new QLabel(" STATUS ");
  statusLabel->setAlignment(Qt::AlignHCenter);
  statusLabel->setMinimumSize(statusLabel->sizeHint());

  messageLabel = new QLabel;
  messageLabel->setIndent(3);

  statusBar()->addWidget(statusLabel);
  statusBar()->addWidget(messageLabel, 1);

//   connect(textBrowser, SIGNAL(currentCellChanged(int, int, int, int)),
// 	  this, SLOT(updateStatusBar()));
//   connect(textBrowser, SIGNAL(modified()),
// 	  this, SLOT(textBrowserModified()));

  updateStatusBar();
}

void mainWindow::updateStatusBar()
{
  statusLabel->setText("");
  messageLabel->setText("");
}

void mainWindow::monitorModified()
{
  setWindowModified(true);
  updateStatusBar();
}

void mainWindow::openDB(){
  int db_id;
  QAction *action = qobject_cast<QAction *>(sender());
  if (action) db_id = action->data().toInt();

  if (db_id < 0){
    statusBar()->showMessage(tr("Opening all connections"), 2000);
    if (pg_dbi::get()->open_connections() != NDB){
      statusBar()->showMessage(tr("Error opening all connections"), 2000);
      return;
    }
    for (int i = 0; i < NDB; ++i) {
      openDBActions[i]->setEnabled(false);
      operateDBActions[i]->setEnabled(true);
      closeDBActions[i]->setEnabled(true);
    }
  } else {
    QString dbName(db_names[db_id]);
    statusBar()->showMessage(tr("Opening connection %1").arg(dbName), 2000);
    if (!pg_dbi::get()->open_connection(db_id)){
      statusBar()->showMessage(tr("Error opening connection to %1").arg(dbName), 2000);
      return;
    }
    openDBActions[db_id]->setEnabled(false);
    operateDBActions[db_id]->setEnabled(true);
    closeDBActions[db_id]->setEnabled(true);
  }

}

void mainWindow::closeDB(){
  int db_id;
  QAction *action = qobject_cast<QAction *>(sender());
  if (action) db_id = action->data().toInt();

  if (db_id < 0){
    statusBar()->showMessage(tr("Closing all connections"), 2000);
    pg_dbi::get()->close_connections();
    for (int i = 0; i < NDB; ++i) {
      openDBActions[i]->setEnabled(true);
      operateDBActions[i]->setEnabled(false);
      closeDBActions[i]->setEnabled(false);
    }
    close();
  } else {
    QString dbName(db_names[db_id]);
    statusBar()->showMessage(tr("Closing connection %1").arg(dbName), 2000);
    pg_dbi::get()->close_connection(db_id);
    openDBActions[db_id]->setEnabled(true);
    operateDBActions[db_id]->setEnabled(false);
    closeDBActions[db_id]->setEnabled(false);
  }

}

bool mainWindow::operateDB(){
  int db_id;
  QAction *action = qobject_cast<QAction *>(sender());
  if (action) db_id = action->data().toInt();

  if (!dbMngDialog) {
    dbMngDialog = new DBDialog(db_id, this);
    for (int i = 0; i < NDB; ++i){
      connect(operateDBActions[i], SIGNAL(triggered()), dbMngDialog, SLOT(setDBId()));
      //connect(dbMngDialog, SIGNAL(), this, SLOT());
    }
  }

  dbMngDialog->show();
  dbMngDialog->raise();
  dbMngDialog->activateWindow();

  return true;
}

void mainWindow::bchMon()
{
  //printf("bad channels monitor triggered\n");
  //openDBActions[2]->activate(QAction::Trigger);
  //if (!pg_dbi::get()->is_opened(2)) return;
  openDBActions[3]->activate(QAction::Trigger);
  if (!pg_dbi::get()->is_opened(3)) return;

  if (!bchMonDialog) {
    bchMonDialog = new BCMDialog(this);
    //connect(, SIGNAL(triggered()), bchMonDialog, SLOT());
  }

  bchMonDialog->show();
  bchMonDialog->raise();
  bchMonDialog->activateWindow();
}

void mainWindow::lmsMon()
{
  //printf("dead in neutrino monitor triggered\n");
  openDBActions[2]->activate(QAction::Trigger);
  if (!pg_dbi::get()->is_opened(2)) return;
  openDBActions[3]->activate(QAction::Trigger);
  if (!pg_dbi::get()->is_opened(3)) return;

  if (!lmsMonDialog) {
    lmsMonDialog = new LMSDialog(this);
    //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
  }

  lmsMonDialog->show();
  lmsMonDialog->raise();
  lmsMonDialog->activateWindow();
}

void mainWindow::mchMon()
{
  //printf("bad channels monitor triggered\n");
  //openDBActions[2]->activate(QAction::Trigger);
  //if (!pg_dbi::get()->is_opened(2)) return;
  openDBActions[5]->activate(QAction::Trigger);
  if (!pg_dbi::get()->is_opened(5)) return;

  if (!mchMonDialog) {
    mchMonDialog = new MCMDialog(this);
    //connect(, SIGNAL(triggered()), bchMonDialog, SLOT());
  }

  mchMonDialog->show();
  mchMonDialog->raise();
  mchMonDialog->activateWindow();
}

void mainWindow::mmsMon()
{
  //printf("dead in neutrino monitor triggered\n");
  openDBActions[2]->activate(QAction::Trigger);
  if (!pg_dbi::get()->is_opened(2)) return;

  if (!mmsMonDialog) {
    mmsMonDialog = new MMSDialog(this);
    //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
  }

  mmsMonDialog->show();
  mmsMonDialog->raise();
  mmsMonDialog->activateWindow();
}

void mainWindow::hvMon()
{
  openAllDBAction->activate(QAction::Trigger);
  for (int i = 0; i < NDB; ++i)
    if (!pg_dbi::get()->is_opened(i)) return;

  if (!hvMonDialog) {
    hvMonDialog = new HVDialog(this);
    //connect(, SIGNAL(triggered()), hvMonDialog, SLOT());
  }

  hvMonDialog->show();
  hvMonDialog->raise();
  hvMonDialog->activateWindow();
}

void mainWindow::dmpMon()
{
  openAllDBAction->activate(QAction::Trigger);
  for (int i = 0; i < NDB; ++i)
    if (!pg_dbi::get()->is_opened(i)) return;

  if (!dmpMonDialog) {
    dmpMonDialog = new DMPDialog(this);
    //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
  }

  dmpMonDialog->show();
  dmpMonDialog->raise();
  dmpMonDialog->activateWindow();
}

void mainWindow::dmmMon()
{
  openAllDBAction->activate(QAction::Trigger);
  for (int i = 0; i < NDB; ++i)
    if (!pg_dbi::get()->is_opened(i)) return;

  if (!dmpmuMonDialog) {
    dmpmuMonDialog = new DMPMuDialog(this);
    //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
  }

  dmpmuMonDialog->show();
  dmpmuMonDialog->raise();
  dmpmuMonDialog->activateWindow();
}

void mainWindow::closeEvent(QCloseEvent *event)
{
//   if (okToContinue()) {
//     writeSettings();
//     event->accept();
//   } else {
//     event->ignore();
//   }
  event->accept();
}

void mainWindow::help()
{
  QMessageBox::about(this, tr("Help"),
		     tr("<h2>\"Help! I need somebody...\"</h2>"
			"<p>                 The Beatles (C)"));
}

void mainWindow::about()
{
  QMessageBox::about(this, tr("dbMon v.0.4.5"),
		     tr("<h2>DataBase MONitor 0.4.5</h2>"
			 "<p>Author(s): Livia Ludhova, Kirill Fomenko"
			 "<p>           for Detector's Stability Group"
		         "<p>Copyright: BOREXINO Collaboration, 2011-2012"));
}

// void mainWindow::chsMon()
// {
//   //printf("dead in neutrino monitor triggered\n");
//   openDBActions[2]->activate(QAction::Trigger);
//   if (!pg_dbi::get()->is_opened(2)) return;
//
//   if (!chsMonDialog) {
//     chsMonDialog = new CHSDialog(this);
//     //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
//   }
//
//   chsMonDialog->show();
//   chsMonDialog->raise();
//   chsMonDialog->activateWindow();
// }
//
// void mainWindow::cbsMon()
// {
//   //printf("dead in neutrino monitor triggered\n");
//   openDBActions[2]->activate(QAction::Trigger);
//   if (!pg_dbi::get()->is_opened(2)) return;
//
//   if (!cbsMonDialog) {
//     cbsMonDialog = new CBSDialog(this);
//     //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
//   }
//
//   cbsMonDialog->show();
//   cbsMonDialog->raise();
//   cbsMonDialog->activateWindow();
// }
//
// void mainWindow::cpsMon()
// {
//   //printf("dead in neutrino monitor triggered\n");
//   openDBActions[2]->activate(QAction::Trigger);
//   if (!pg_dbi::get()->is_opened(2)) return;
//
//   if (!cpsMonDialog) {
//     cpsMonDialog = new CPSDialog(this);
//     //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
//   }
//
//   cpsMonDialog->show();
//   cpsMonDialog->raise();
//   cpsMonDialog->activateWindow();
// }
//
// void mainWindow::tmsMon()
// {
//   //printf("dead in neutrino monitor triggered\n");
//   openDBActions[2]->activate(QAction::Trigger);
//   if (!pg_dbi::get()->is_opened(2)) return;
//
//   if (!tmsMonDialog) {
//     tmsMonDialog = new TMSDialog(this);
//     //connect(, SIGNAL(triggered()), lmsMonDialog, SLOT());
//   }
//
//   tmsMonDialog->show();
//   tmsMonDialog->raise();
//   tmsMonDialog->activateWindow();
// }
