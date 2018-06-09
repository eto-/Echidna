#include <QtGui>

#include "mainWindow.h"
#include "mchmondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QButtonGroup>
#include <QRadioButton>
#include <QPushButton>
#include <QTableWidget>
#include <QCheckBox>
#include <QStringList>

// *********************** outline window

MCMDialog::MCMDialog(QWidget *parent)
  : QDialog(parent), run_n(0), tt1_rt(0), tt2_rt(0), mean_rt(0)  //, bad_ls(0), off_ls(0)
{
  str = new QString();

  //bad_ls = new int[CH_NUM]; off_ls = new int[CH_NUM];

  //QRegExp regExp("[A-Za-z][1-9][0-9]{1,3}");
  QRegExp regExp("[0-9]{0,4}");
  stepLabel = new QLabel(tr("&step value"), this);
  stepEdit  = new QLineEdit(this);
  stepLabel->setBuddy(stepEdit);
  stepEdit->setText(tr("100"));
  stepEdit->setValidator(new QRegExpValidator(regExp, this));

  regExp.setPattern("[0-9]{0,6}");
  run1Label = new QLabel(tr("&first Run"), this);
  run1Edit  = new QLineEdit(this);
  run1Label->setBuddy(run1Edit);
    //run1Edit->setText(tr("7733"));
  run1Edit->setText(tr("%1").arg(mainWindow::current_run()-700));
  run1Edit->setValidator(new QRegExpValidator(regExp, this));

  run2Label = new QLabel(tr("&last Run"), this);
  run2Edit  = new QLineEdit(this);
  run2Label->setBuddy(run2Edit);
  run2Edit->setText(tr("%1").arg(mainWindow::current_run()));
  run2Edit->setValidator(new QRegExpValidator(regExp, this));

  connect(stepEdit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));
  connect(run1Edit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));
  connect(run2Edit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));

  prevRadio = new QRadioButton(tr("to &previous           "));
  averRadio = new QRadioButton(tr("to &average            "));
  prevRadio->setChecked(true);
  QButtonGroup *methRadGrp = new QButtonGroup(this);
  methRadGrp->addButton(prevRadio, 1);
  methRadGrp->addButton(averRadio, 2);

  tt1Radio = new QRadioButton(tr("TT&1"));
  tt2Radio = new QRadioButton(tr("TT&2"));
  tt1Radio->setChecked(true);
  QButtonGroup *typeRadGrp = new QButtonGroup(this);
  typeRadGrp->addButton(tt1Radio, 1);
  typeRadGrp->addButton(tt2Radio, 2);

  regExp.setPattern("[0-9]{1,4}");
  meanRadio = new QCheckBox(tr("&mean over"));
  meanRadio->setChecked(true);
  meanEdit  = new QLineEdit(this);
  meanEdit->setValidator(new QRegExpValidator(regExp, this));
  meanEdit->setText(tr("50"));
  meanEdit->setFixedWidth((int)(0.4*meanEdit->sizeHint().width()));
  meanEdit->setFixedHeight((int)(1.0*meanEdit->sizeHint().height()));
  meanLabel = new QLabel(tr("        "), this);
  meanLabel->setBuddy(meanEdit);

  //connect(meanRadio, SIGNAL(toggled(bool)), meanEdit,  SLOT(setEnabled(bool)));
  //connect(meanRadio, SIGNAL(toggled(bool)), meanLabel, SLOT(setEnabled(bool)));

  proceedButton = new QPushButton(tr("P&roceed"), this);
  proceedButton->setEnabled(false);

  cancelButton = new QPushButton(tr("Can&cel"), this);
  cancelButton->setDefault(true);
  cancelButton->setEnabled(true);

  connect(proceedButton, SIGNAL(clicked()), this, SLOT(showChannels()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));

  queryContent = new QTableWidget(10, 4, this);
  QStringList tableLabels;
  tableLabels.append(tr(" Run number "));
  tableLabels.append(tr("   TT1 rate   "));
  tableLabels.append(tr("   TT2 rate   "));
  tableLabels.append(tr("mean over group"));
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
  paramsLayout->addLayout(run1Layout);
  paramsLayout->insertStretch(2, 1);
  paramsLayout->addLayout(stepLayout);
  paramsLayout->insertStretch(4, 1);
  paramsLayout->addLayout(run2Layout);
  //paramsLayout->insertStretch(4, 1);

  QHBoxLayout *radioButtonsLayout11 = new QHBoxLayout;
  radioButtonsLayout11->addWidget(prevRadio);
  QHBoxLayout *radioButtonsLayout12 = new QHBoxLayout;
  radioButtonsLayout12->addWidget(averRadio);
  QVBoxLayout *radioButtonsLayout1 = new QVBoxLayout;
  radioButtonsLayout1->addLayout(radioButtonsLayout11);
  radioButtonsLayout1->addLayout(radioButtonsLayout12);

  QHBoxLayout *radioButtonsLayout21 = new QHBoxLayout;
  radioButtonsLayout21->addWidget(tt1Radio, 0);
  radioButtonsLayout21->addWidget(tt2Radio, 0);
  radioButtonsLayout21->insertStretch(3, 1);
  QHBoxLayout *radioButtonsLayout22 = new QHBoxLayout;
  radioButtonsLayout22->addWidget(meanRadio, 0);
  radioButtonsLayout22->addWidget(meanEdit, 0);
  radioButtonsLayout22->addWidget(meanLabel, 0);
  radioButtonsLayout22->insertStretch(4, 1);
  QVBoxLayout *radioButtonsLayout2 = new QVBoxLayout;
  radioButtonsLayout2->addLayout(radioButtonsLayout21);
  radioButtonsLayout2->addLayout(radioButtonsLayout22);

  QVBoxLayout *controlButtonsLayout = new QVBoxLayout;
  controlButtonsLayout->addWidget(proceedButton);
  controlButtonsLayout->addWidget(cancelButton);

  QHBoxLayout *buttonsLayout = new QHBoxLayout;
  buttonsLayout->addLayout(radioButtonsLayout2);
  //buttonsLayout->insertStretch(2, 1);
  buttonsLayout->addLayout(radioButtonsLayout1);
  buttonsLayout->insertStretch(3, 1);
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
  setWindowTitle(tr("Muon Trigger rates"));
}

MCMDialog::~MCMDialog(){
  if (run_n)  delete[] run_n;
  if (tt1_rt) delete[] tt1_rt;
  if (tt2_rt) delete[] tt2_rt;
  if (mean_rt) delete[] mean_rt;
//   if (bad_ls) delete[] bad_ls;
//   if (off_ls) delete[] off_ls;
}

void MCMDialog::enableProceedButton()
{
  if (!stepEdit->text().isEmpty() && !run1Edit->text().isEmpty() && !run2Edit->text().isEmpty() && !meanEdit->text().isEmpty()){
    proceedButton->setEnabled(true);
    cancelButton->setDefault(false);
    proceedButton->setDefault(true);
  } else{
    proceedButton->setEnabled(false);
    proceedButton->setDefault(false);
    cancelButton->setDefault(true);
  }
}

void MCMDialog::showChannels()
{
  run1 = run1Edit->text().toInt();
  run2 = run2Edit->text().toInt();
  step = stepEdit->text().toInt();
  grpn = meanEdit->text().toInt();

  run2 = run2 > run1 ? run2 : run1;
  run2Edit->setText(tr("%1").arg(run2));

  int tt1 = tt1Radio->isChecked()  ?  1 : 0;
  int grp = meanRadio->isChecked() ?  1 : 0;

  if (!run1 || !run2 || !step || (grp && !grpn)) return;

  sprintf(the_query, "SELECT \"RunNumber\", \"MtbTt1Rate\", \"MtbTt2Rate\" FROM \"ValidRuns\" \
                      WHERE \"RunNumber\" > %d AND \"RunNumber\" < %d ORDER BY \"RunNumber\" ", run1-1, run2+1);
  pg_dbi::get()->query(BX_RUNVALIDATION, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  if (run_n)  delete[] run_n;
  run_n  = new int[lines];
  if (tt1_rt) delete[] tt1_rt;
  tt1_rt = new float[lines];
  if (tt2_rt) delete[] tt2_rt;
  tt2_rt = new float[lines];
  if (mean_rt) delete[] mean_rt;
  mean_rt = new float[lines];

  int lines_to_show = 0, groups = 0;
  float aver_value = 0., prev_value = 0., mean_over_group = 0.;

  for (int i = 0; i < lines; ++i){
    if (!QString::fromAscii(pg_dbi::get()->query_value(i, 2)).toInt()) continue; // skip empty lines
    run_n[lines_to_show]   = QString::fromAscii(pg_dbi::get()->query_value(i, 0)).toInt();
    tt1_rt[lines_to_show]  = QString::fromAscii(pg_dbi::get()->query_value(i, 1)).toFloat();
    tt2_rt[lines_to_show]  = QString::fromAscii(pg_dbi::get()->query_value(i, 2)).toFloat();
    mean_rt[lines_to_show] = 0.;
    if (lines_to_show && !(lines_to_show % grpn)){
      mean_rt[lines_to_show] = mean_over_group/grpn;
      if (grp) aver_value += mean_over_group/grpn;
      mean_over_group = 0.; ++groups;
    }
    if (tt1){
      mean_over_group += tt1_rt[lines_to_show];
      if (!grp) aver_value += tt1_rt[lines_to_show];
    }
    else {
      mean_over_group += tt2_rt[lines_to_show];
      if (!grp) aver_value += tt2_rt[lines_to_show];
    }
    ++lines_to_show;
  }

  if (lines_to_show){
    if (grp){
      if (groups){
        aver_value /= groups;
        prev_value = mean_rt[grpn];
      }
    } else {
      aver_value /= lines_to_show;
      if (tt1) prev_value = tt1_rt[0];
      else     prev_value = tt2_rt[0];
    }
  }

  queryContent->clearContents();
  if (grp) queryContent->setRowCount(groups);
  else queryContent->setRowCount(lines_to_show);

  int table_line = -1;
  for (int i = 0; i < lines_to_show; ++i){
    if (grp){
      if (!i || (i % grpn)) continue;
      ++table_line;
    } else
      table_line = i;

    QTableWidgetItem *tableItems0 = new QTableWidgetItem(QString::number(run_n[i]));
    tableItems0->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems0->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(table_line, 0, tableItems0);

    QTableWidgetItem *tableItems1 = new QTableWidgetItem(QString::number(tt1_rt[i]));
    tableItems1->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(table_line, 1, tableItems1);

    QTableWidgetItem *tableItems2 = new QTableWidgetItem(QString::number(tt2_rt[i]));
    tableItems2->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(table_line, 2, tableItems2);

    QTableWidgetItem *tableItems3 = new QTableWidgetItem(QString::number(mean_rt[i]));
    tableItems3->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems3->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(table_line, 3, tableItems3);
  }

  mainWindow::get()->dispatchMsg(str->sprintf(" >>> Muon trigger rates: runs %d - %d\n", run1, run2), 0);
  mainWindow::get()->dispatchMsg(str->sprintf(" step: %d, method: ", step), 0);
  if (prevRadio->isChecked())
    mainWindow::get()->dispatchMsg(str->sprintf("with respect to previous, reference value: "), 0);
  else
    mainWindow::get()->dispatchMsg(str->sprintf("with respect to average, reference value: "), 0);
  if (grp){
    if (tt1)
      mainWindow::get()->dispatchMsg(str->sprintf("TT1 rate, mean over %d runs\n", grpn), 0);
    else
      mainWindow::get()->dispatchMsg(str->sprintf("TT2 rate, mean over %d runs\n", grpn), 0);
    mainWindow::get()->dispatchMsg(str->sprintf("     run     TT1     TT2    mean\n"), 0);
  } else {
    if (tt1)
      mainWindow::get()->dispatchMsg(str->sprintf("TT1 rate\n"), 0);
    else
      mainWindow::get()->dispatchMsg(str->sprintf("TT2 rate\n"), 0);
    mainWindow::get()->dispatchMsg(str->sprintf("     run     TT1     TT2\n"), 0);
  }

  if (!grp) grpn = 1;

  QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);
  int match_found, total_found = 0;
  table_line = -1;
  for (int i = 0; i < lines_to_show; i += grpn){
    match_found = 0;
    if (grp){
      if (!i) continue;
      if ( (averRadio->isChecked() && (mean_rt[i] > (aver_value + step) || mean_rt[i] < (aver_value - step))) ||
	   (prevRadio->isChecked() && (mean_rt[i] > (prev_value + step) || mean_rt[i] < (prev_value - step))) ) match_found = 1;
      prev_value = mean_rt[i];
      ++table_line;
    } else {
      if (tt1){
	if ( (averRadio->isChecked() && (tt1_rt[i] > (aver_value + step) || tt1_rt[i] < (aver_value - step))) ||
	     (prevRadio->isChecked() && (tt1_rt[i] > (prev_value + step) || tt1_rt[i] < (prev_value - step))) ) match_found = 1;
	prev_value = tt1_rt[i];
      } else {
	if ( (averRadio->isChecked() && (tt2_rt[i] > (aver_value + step) || tt2_rt[i] < (aver_value - step))) ||
	     (prevRadio->isChecked() && (tt2_rt[i] > (prev_value + step) || tt2_rt[i] < (prev_value - step))) ) match_found = 1;
	prev_value = tt2_rt[i];
      }
      table_line = i;
    }
    if (match_found){
      queryContent->item(table_line, 0)->setBackground(bg_col); queryContent->item(table_line, 0)->setForeground(fg_col);
      queryContent->item(table_line, 1)->setBackground(bg_col); queryContent->item(table_line, 1)->setForeground(fg_col);
      queryContent->item(table_line, 2)->setBackground(bg_col); queryContent->item(table_line, 2)->setForeground(fg_col);
      queryContent->item(table_line, 3)->setBackground(bg_col); queryContent->item(table_line, 3)->setForeground(fg_col);

      // printing
      if (grp)
	mainWindow::get()->dispatchMsg(str->sprintf(" %7d %7.0f %7.0f %7.0f\n", run_n[i], tt1_rt[i], tt2_rt[i], mean_rt[i]), 0);
      else
	mainWindow::get()->dispatchMsg(str->sprintf(" %7d %7.0f %7.0f\n", run_n[i], tt1_rt[i], tt2_rt[i]), 0);
      ++total_found;
    }
  }
  mainWindow::get()->dispatchMsg(str->sprintf(" total: %d", total_found), 0);
  if (grp && averRadio->isChecked()) mainWindow::get()->dispatchMsg(str->sprintf("  (average value over groups: %.2f)\n", aver_value), 0);
  else mainWindow::get()->dispatchMsg(str->sprintf("\n"), 0);

  QStringList tableLabels;
  for (int i = 0; i < lines_to_show; i += grpn) tableLabels.append(tr("#"));
  queryContent->setVerticalHeaderLabels(tableLabels);

  queryContent->resizeColumnsToContents();
  //queryContent->adjustSize();
  queryContent->update();
  queryContent->show();
  //update();

  fflush(stdout);
}

