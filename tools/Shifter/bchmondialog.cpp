#include <QtGui>

#include "mainWindow.h"
#include "bchmondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QRadioButton>
#include <QPushButton>
#include <QTableWidget>
#include <QStringList>

// *********************** outline window

BCMDialog::BCMDialog(QWidget *parent)
  : QDialog(parent), run_n(0), bad_ch(0), off_ch(0), sum_ch(0), bad_ls(0), off_ls(0)
{
  str = new QString();

  bad_ls = new int[CH_NUM]; off_ls = new int[CH_NUM];

  //QRegExp regExp("[A-Za-z][1-9][0-9]{1,3}");
  QRegExp regExp("[0-9]{0,4}");

  stepLabel = new QLabel(tr("value &step"), this);
  stepEdit  = new QLineEdit(this);
  stepLabel->setBuddy(stepEdit);
  stepEdit->setText(tr("5"));
  stepEdit->setValidator(new QRegExpValidator(regExp, this));

  regExp.setPattern("[0-9]{0,6}");
  run1Label = new QLabel(tr("&first Run"), this);
  run1Edit  = new QLineEdit(this);
  run1Label->setBuddy(run1Edit);
    //run1Edit->setText(tr("7733"));
  run1Edit->setText(tr("%1").arg(mainWindow::current_run()-100));
  run1Edit->setValidator(new QRegExpValidator(regExp, this));

  run2Label = new QLabel(tr("&last Run"), this);
  run2Edit  = new QLineEdit(this);
  run2Label->setBuddy(run2Edit);
  run2Edit->setText(tr("%1").arg(mainWindow::current_run()));
  run2Edit->setValidator(new QRegExpValidator(regExp, this));

  connect(stepEdit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));
  connect(run1Edit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));
  connect(run2Edit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));

  prevRadio = new QRadioButton(tr("to &previous"), this);
  averRadio = new QRadioButton(tr("to &average "), this);
  prevRadio->setChecked(true);

  regExp.setPattern("[0-9]{1,2}");
  pcntLabel = new QLabel(tr(" &threshold, %"), this);
  pcntEdit  = new QLineEdit(this);
  pcntLabel->setBuddy(pcntEdit);
  pcntEdit->setText(tr("80"));
  pcntEdit->setValidator(new QRegExpValidator(regExp, this));
  pcntEdit->setFixedWidth((int)(0.3*pcntEdit->sizeHint().width()));

  proceedButton = new QPushButton(tr("P&roceed"), this);
  proceedButton->setEnabled(false);

  cancelButton = new QPushButton(tr("Can&cel"), this);
  cancelButton->setDefault(true);
  cancelButton->setEnabled(true);

  connect(proceedButton, SIGNAL(clicked()), this, SLOT(showChannels()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));

  queryContent = new QTableWidget(10, 4, this);
  QStringList tableLabels;
  tableLabels.append(tr("Run number"));
  tableLabels.append(tr("Bad channels"));
  tableLabels.append(tr("Off channels"));
  tableLabels.append(tr("  Total  "));
  queryContent->setHorizontalHeaderLabels(tableLabels);
  queryContent->resizeColumnsToContents();

  QVBoxLayout *run1Layout = new QVBoxLayout;
  run1Layout->addWidget(run1Label);
  run1Layout->addWidget(run1Edit);

  QVBoxLayout *run2Layout = new QVBoxLayout;
  run2Layout->addWidget(run2Label);
  run2Layout->addWidget(run2Edit);

  QVBoxLayout *stepLayout = new QVBoxLayout;
  stepLayout->addWidget(stepLabel);
  stepLayout->addWidget(stepEdit);

  QHBoxLayout *paramsLayout = new QHBoxLayout;
  paramsLayout->addLayout(stepLayout);
  paramsLayout->insertStretch(2, 1);
  paramsLayout->addLayout(run1Layout);
  paramsLayout->insertStretch(4, 1);
  paramsLayout->addLayout(run2Layout);
  //paramsLayout->insertStretch(4, 1);

  QVBoxLayout *radioButtonsLayout = new QVBoxLayout;
  radioButtonsLayout->addWidget(prevRadio);
  radioButtonsLayout->addWidget(averRadio);

  QHBoxLayout *thhLayout = new QHBoxLayout;
  thhLayout->addWidget(pcntEdit);
  thhLayout->addWidget(pcntLabel);

//   QVBoxLayout *thresholdLayout = new QVBoxLayout;
//   thresholdLayout->addLayout(thhLayout);
//   thresholdLayout->insertStretch(2,1);

  QVBoxLayout *controlButtonsLayout = new QVBoxLayout;
  controlButtonsLayout->addWidget(proceedButton);
  controlButtonsLayout->addWidget(cancelButton);

  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  buttonsLayout->addLayout(radioButtonsLayout);
  buttonsLayout->insertStretch(2, 1);
  buttonsLayout->addLayout(thhLayout);
  buttonsLayout->insertStretch(4, 1);
  buttonsLayout->addLayout(controlButtonsLayout);

  QVBoxLayout *tableLayout = new QVBoxLayout;
  tableLayout->addWidget(queryContent, 1);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(paramsLayout, 0);
  mainLayout->addLayout(buttonsLayout, 0);
  mainLayout->addLayout(tableLayout, 1);
  //mainLayout->insertStretch(4, 1);
  setLayout(mainLayout);

  //paramsLayout->setFixedHeight(paramsLayout->sizeHint().height());
  //buttonsLayout->setFixedHeight(bottonsLayout->sizeHint().height());

  setMinimumWidth((int)(1.05*sizeHint().width()));
  setMaximumWidth((int)(1.1*sizeHint().width()));
  //setFixedWidth((int)(1.05*sizeHint().width()));

  queryContent->hide();
  setWindowTitle(tr("Bad/Off Channels stats"));
}

BCMDialog::~BCMDialog(){
  if (run_n)  delete[] run_n;
  if (bad_ch) delete[] bad_ch;
  if (off_ch) delete[] off_ch;
  if (sum_ch) delete[] sum_ch;
  if (bad_ls) delete[] bad_ls;
  if (off_ls) delete[] off_ls;
}

void BCMDialog::enableProceedButton()
{
  if (!stepEdit->text().isEmpty() && !run1Edit->text().isEmpty() && !run2Edit->text().isEmpty() && !pcntEdit->text().isEmpty()){
    proceedButton->setEnabled(true);
    cancelButton->setDefault(false);
    proceedButton->setDefault(true);
  } else{
    proceedButton->setEnabled(false);
    proceedButton->setDefault(false);
    cancelButton->setDefault(true);
  }
}

void BCMDialog::showChannels()
{

  float aver_value = 0;
  int prev_value = 0;

  run1 = run1Edit->text().toInt();
  run2 = run2Edit->text().toInt();
  step = stepEdit->text().toInt();
  pcnt = pcntEdit->text().toInt();

  run2 = run2 > run1 ? run2 : run1;
  run2Edit->setText(tr("%1").arg(run2));

  if (!run1 || !run2 || !step || !pcnt) return;

  sprintf(the_query, "SELECT \"RunNumber\", \"BadChannelsNumber\", \"OffChannelsNumber\", \"BadChannelsList\", \"OffChannelsList\" FROM \"LabenPrecalibDecodingQuality\" \
                      WHERE \"RunNumber\" > %d AND \"RunNumber\" < %d ORDER BY \"RunNumber\" ", run1-1, run2+1);
  pg_dbi::get()->query(BX_PRECALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  if (run_n)  delete[] run_n;
  run_n  = new int[lines];
  if (bad_ch) delete[] bad_ch;
  bad_ch = new int[lines];
  if (off_ch) delete[] off_ch;
  off_ch = new int[lines];
  if (sum_ch) delete[] sum_ch;
  sum_ch = new int[lines];
  for (int i = 0; i <= CH_NUM; ++i){
    bad_ls[i] = off_ls[i] = 0;
  }

  //return;
  int ch;
  QRegExp regExpParser(",");
  QRegExp regExpPattern("[{}]");
  QStringList bad_ch_list, off_ch_list;
  for (int i = 0; i < lines; ++i){
    run_n[i]  = QString::fromAscii(pg_dbi::get()->query_value(i, 0)).toInt();
    bad_ch[i] = QString::fromAscii(pg_dbi::get()->query_value(i, 1)).toInt();
    off_ch[i] = QString::fromAscii(pg_dbi::get()->query_value(i, 2)).toInt();
    sum_ch[i] = bad_ch[i] + off_ch[i];
    aver_value += sum_ch[i];
    bad_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(i, 3)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);
    off_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(i, 4)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);
    //bad_ch_list.removeFirst(); bad_ch_list.removeLast(); off_ch_list.removeFirst(); off_ch_list.removeLast();
    //printf("Run %d bad channels:", run_n[i]);
    for (int j = 0; j < bad_ch_list.size(); ++j){
      ch = bad_ch_list.at(j).toInt();
      if (ch < CH_NUM+1 && ch > 0) ++bad_ls[ch-1];
    }
    for (int j = 0; j < off_ch_list.size(); ++j){
      ch = off_ch_list.at(j).toInt();
      if (ch < CH_NUM+1 && ch > 0) ++off_ls[ch-1];
    }
    //printf(" %s", bad_ch_list.at(j).toLocal8Bit().constData());
    bad_ch_list.clear(); off_ch_list.clear();
    //for (int j = 1; j <= CH_NUM; ++j){
      //if (QString::fromAscii(pg_dbi::get()->query_value(i, 3)).contains(QString::number(j), Qt::CaseInsensitive)) ++bad_ls[j];
      //if (QString::fromAscii(pg_dbi::get()->query_value(i, 4)).contains(QString::number(j), Qt::CaseInsensitive)) ++off_ls[j];
    //}
  }
  if (lines){
    prev_value = sum_ch[0]; aver_value /= lines;
  }

  queryContent->clearContents();
  queryContent->setRowCount(lines);

  for (int i = 0; i < lines; ++i){
    QTableWidgetItem *tableItems0 = new QTableWidgetItem(QString::number(run_n[i]));
    tableItems0->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems0->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(i, 0, tableItems0);

    QTableWidgetItem *tableItems1 = new QTableWidgetItem(QString::number(bad_ch[i]));
    tableItems1->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(i, 1, tableItems1);

    QTableWidgetItem *tableItems2 = new QTableWidgetItem(QString::number(off_ch[i]));
    tableItems2->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(i, 2, tableItems2);

    QTableWidgetItem *tableItems3 = new QTableWidgetItem(QString::number(sum_ch[i]));
    tableItems3->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems3->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(i, 3, tableItems3);
  }

  mainWindow::get()->dispatchMsg(str->sprintf(" >>> Bad/Off channels in precalibration: runs %d - %d\n", run1, run2), 0);
  mainWindow::get()->dispatchMsg(str->sprintf(" step: %d, method: ", step), 0);
  if (prevRadio->isChecked())
    mainWindow::get()->dispatchMsg(str->sprintf("with respect to previous\n"), 0);
  else
    mainWindow::get()->dispatchMsg(str->sprintf("with respect to average\n"), 0);
  mainWindow::get()->dispatchMsg(str->sprintf("     run   N_bad   N_off   total\n"), 0);

  QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);
  int total_found = 0;
  for (int i = 0; i < lines; ++i){
    if ( (averRadio->isChecked() && sum_ch[i] >= (aver_value + step)) || (prevRadio->isChecked() && sum_ch[i] >= (prev_value + step)) ){
      queryContent->item(i, 0)->setBackground(bg_col); queryContent->item(i, 0)->setForeground(fg_col);
      queryContent->item(i, 1)->setBackground(bg_col); queryContent->item(i, 1)->setForeground(fg_col);
      queryContent->item(i, 2)->setBackground(bg_col); queryContent->item(i, 2)->setForeground(fg_col);
      queryContent->item(i, 3)->setBackground(bg_col); queryContent->item(i, 3)->setForeground(fg_col);

      // printing
      mainWindow::get()->dispatchMsg(str->sprintf(" %7d %7d %7d %7d\n", run_n[i], bad_ch[i], off_ch[i], sum_ch[i]), 0);
      ++total_found;
    }
    prev_value = sum_ch[i];
  }
  mainWindow::get()->dispatchMsg(str->sprintf(" total: %d\n", total_found), 0);

  QStringList tableLabels;
  for (int i = 0; i < lines; ++i) tableLabels.append(tr("#"));
  queryContent->setVerticalHeaderLabels(tableLabels);

  float bad_ocr, off_ocr;
  mainWindow::get()->dispatchMsg(str->sprintf(" failed praclib in more then %d%% of selected runs (%d - %d):\n", pcnt, run1, run2), 0);
  for (int i = 0; i < CH_NUM; ++i){
    bad_ocr = 100.*((float)bad_ls[i])/((float)lines);
    off_ocr = 100.*((float)off_ls[i])/((float)lines);
    if ((bad_ocr + off_ocr) > pcnt)
      mainWindow::get()->dispatchMsg(str->sprintf("   lg %4d: %3.0f%% of runs, bad: %3.0f%%, off: %3.0f%%\n", i+1, bad_ocr + off_ocr, bad_ocr, off_ocr), 0);
  }

  queryContent->resizeColumnsToContents();
    //queryContent->adjustSize();
  queryContent->update();
  queryContent->show();
  //update();

  //fflush(stdout);
}

/*
  // disconnected channels
  int dsc_lines, dsc_columns;
  char ch_state[1024];
  for (int i = 0; i < lines; ++i){
    sprintf(the_query, "SELECT \"ChannelID\",\"RunNumber\",\"Disconnected\" FROM \"DisconnectedPmts\" \
	                WHERE \"RunNumber\"<%d ORDER BY \"ChannelID\",\"RunNumber\" ", run_n[i]+1);
    pg_dbi::get()->test_query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(dsc_lines, dsc_columns);

    for (int j = 0; j < dsc_lines; ++j){
      ch = atoi(pg_dbi::get()->query_value(j, 0));
      sprintf(ch_state, pg_dbi::get()->query_value(j, 2));
      if (!strcmp(ch_state, "t")){
	--off_ch[i];
	--sum_ch[i];
	--off_ls[ch-1];
      }
    }
  }  */
