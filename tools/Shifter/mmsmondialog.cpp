#include <QtGui>

#include "mainWindow.h"
#include "mmsmondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QTableWidget>
#include <QStringList>
#include <QCheckBox>

const char MMSDialog::logic[][1024] = {"OR", "AND", "AND NOT"};

const char MMSDialog::fields[][1024] =
{"", "Multiplicity"}; //, "ChargeStatus", "ChargeBaseStatus", "ChargePeakStatus", "TimingStatus"};

const char MMSDialog::types1[][1024] =
{ "",
  "dead_in_muon",         "dead_in_neutrino",         "dead_in_pulser",         "dead_in_laser",
  "low_eff_in_muon",      "low_eff_in_neutrino",      "low_eff_in_pulser",      "low_eff_in_laser",
  "hot_in_muon",          "hot_in_neutrino",          "hot_in_pulser",          "hot_in_laser",
  "retriggering_in_muon", "retriggering_in_neutrino", "retriggering_in_pulser", "retriggering_in_laser"
};

// const char MMSDialog::types1[][1024] =
// { "",
//   "dead_in_neutrino", "dead_in_pulser", "dead_in_laser", "dead_in_raw",
//   "low_eff_in_neutrino", "low_eff_in_pulser", "low_eff_in_laser", "low_eff_in_raw",
//   "hot_in_neutrino", "hot_in_pulser", "hot_in_laser", "hot_in_raw",
//   "retriggering_in_neutrino",  "retriggering_in_pulser", "retriggering_in_laser",
//   "loosing_raw_hits_in_neutrino", "loosing_raw_hits_in_pulser", "loosing_raw_hits_in_laser",
//   "fifo_empty", "fifo_full"
// };

const char MMSDialog::types2[][1024] =
{ "", "many_zero", "many_FF", "too_spread", "shifted_from_mean", "bad_rms", "very_small_rms", "many_negative_values", "low_gain", "high_gain" };

const char MMSDialog::types3[][1024] =
{ "", "laser_ref_correl", "trigger_ref_correl_in_pulser", "trigger_ref_correl_in_laser", "trigger_ref_correl_in_neutrino",
  "end_of_gate_correl_in_pulser", "end_of_gate_correl_in_laser", "end_of_gate_correl_in_neutrino", "bad_timing_shape_in_laser"
};

// *********************** outline window

MMSDialog::MMSDialog(QWidget *parent)
  : QDialog(parent), field1(0), field2(0), field3(0), type1(0), type2(0), type3(0)
{
  str = new QString();

  for (int i = 0; i < CH_NUM; ++i){
    type1_values[i] = new char[4096]; type2_values[i] = new char[4096]; type3_values[i] = new char[4096];
  }
  dsc_ch      = new int[CH_NUM];
  peak_values = new float[CH_NUM];
  mean_values = new float[CH_NUM];
  sig_values  = new float[CH_NUM];
  rms_values  = new float[CH_NUM];

  QRegExp regExp("[0-9]{0,6}");
  //regExp.setPattern("[1-9]{1,1}[0-9]{0,4}");
  runLabel = new QLabel(tr("RUN &Number:"), this);
  runEdit  = new QLineEdit(this);
  runLabel->setBuddy(runEdit);
  //runEdit->setText(tr("7733"));
  runEdit->setText(tr("%1").arg(mainWindow::current_run()));
  runEdit->setValidator(new QRegExpValidator(regExp, this));

  connect(runEdit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));

  selectionBoxFld1 = new QComboBox(this);
  for (int i = 0; i < N_fields; ++i) selectionBoxFld1->insertItem(i, QString::fromAscii(fields[i]));
  selectionBoxTyp1 = new QComboBox(this);
  //selectionBoxTyp1->setSizeAdjustPolicy(QComboBox::AdjustToContents);

  connect(selectionBoxFld1, SIGNAL(currentIndexChanged(int)), this, SLOT(setField1(int)));
  connect(selectionBoxTyp1, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(setType1(const QString &)));

  selectionBoxFld2 = new QComboBox(this);
  for (int i = 0; i < N_fields; ++i) selectionBoxFld2->insertItem(i, QString::fromAscii(fields[i]));
  selectionBoxTyp2 = new QComboBox(this);
  //selectionBoxTyp2->setSizeAdjustPolicy(QComboBox::AdjustToContents);

  connect(selectionBoxFld2, SIGNAL(currentIndexChanged(int)), this, SLOT(setField2(int)));
  connect(selectionBoxTyp2, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(setType2(const QString &)));

  selectionBoxFld3 = new QComboBox(this);
  for (int i = 0; i < N_fields; ++i) selectionBoxFld3->insertItem(i, QString::fromAscii(fields[i]));
  selectionBoxTyp3 = new QComboBox(this);
  //selectionBoxTyp3->setSizeAdjustPolicy(QComboBox::AdjustToContents);

  connect(selectionBoxFld3, SIGNAL(currentIndexChanged(int)), this, SLOT(setField3(int)));
  connect(selectionBoxTyp3, SIGNAL(currentIndexChanged(const QString &)), this, SLOT(setType3(const QString &)));

  selectionBoxLgc1 = new QComboBox(this);
  for (int i = 0; i < N_conds; ++i) selectionBoxLgc1->insertItem(i, QString::fromAscii(logic[i]));

  selectionBoxLgc2 = new QComboBox(this);
  for (int i = 0; i < N_conds; ++i) selectionBoxLgc2->insertItem(i, QString::fromAscii(logic[i]));

  proceedButton = new QPushButton(tr("P&roceed"), this);
  proceedButton->setEnabled(false);

  cancelButton = new QPushButton(tr("Can&cel"), this);
  cancelButton->setDefault(true);
  cancelButton->setEnabled(true);

  connect(proceedButton, SIGNAL(clicked()), this, SLOT(showChannels()));
  connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));

  queryContent = new QTableWidget(0, 4, this);
  queryContent->hide();

  boxLabel = new QLabel(tr("show all (bad/off/disconnected excluded)"), this);
  checkBox = new QCheckBox(this);
  checkBox->setChecked(false);

  prnLabel = new QLabel(tr("print out disconnected"), this);
  prnCheckBox = new QCheckBox(this);
  prnCheckBox->setChecked(false);

  enableProceedButton();

  QVBoxLayout *runLayout = new QVBoxLayout;
  runLayout->addWidget(runLabel);
  runLayout->addWidget(runEdit);

  QVBoxLayout *comboLayout1 = new QVBoxLayout;
  comboLayout1->addWidget(selectionBoxFld1);
  //comboLayout->insertStretch(1, 1);
  comboLayout1->addWidget(selectionBoxTyp1);

  QVBoxLayout *comboLayout11 = new QVBoxLayout;
  comboLayout11->addWidget(selectionBoxLgc1);

  QVBoxLayout *comboLayout2 = new QVBoxLayout;
  comboLayout2->addWidget(selectionBoxFld2);
  //comboLayout->insertStretch(1, 1);
  comboLayout2->addWidget(selectionBoxTyp2);

  QVBoxLayout *comboLayout22 = new QVBoxLayout;
  comboLayout22->addWidget(selectionBoxLgc2);

  QVBoxLayout *comboLayout3 = new QVBoxLayout;
  comboLayout3->addWidget(selectionBoxFld3);
  //comboLayout->insertStretch(1, 1);
  comboLayout3->addWidget(selectionBoxTyp3);

  QVBoxLayout *buttonsLayout = new QVBoxLayout;
  buttonsLayout->addWidget(proceedButton);
  buttonsLayout->addWidget(cancelButton);

  QHBoxLayout *paramsLayout = new QHBoxLayout;
  paramsLayout->addLayout(runLayout);
  paramsLayout->addLayout(comboLayout1);
  paramsLayout->addLayout(comboLayout11);
  paramsLayout->addLayout(comboLayout2);
  paramsLayout->addLayout(comboLayout22);
  paramsLayout->addLayout(comboLayout3);
  paramsLayout->addLayout(buttonsLayout);
  paramsLayout->insertStretch(6, 1);

  QVBoxLayout *tableLayout = new QVBoxLayout;
  tableLayout->addWidget(queryContent, 1);

  QHBoxLayout *boxLayout = new QHBoxLayout;
  boxLayout->addWidget(checkBox);
  boxLayout->addWidget(boxLabel);
  boxLayout->insertStretch(3, 1);
  boxLayout->addWidget(prnLabel);
  boxLayout->addWidget(prnCheckBox);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(paramsLayout);
  mainLayout->addLayout(tableLayout, 1);
  mainLayout->addLayout(boxLayout);
  setLayout(mainLayout);

  setMinimumWidth(paramsLayout->sizeHint().width());
  //setFixedWidth(sizeHint().width());

  setWindowTitle(tr("Muon Channels status"));

  preInit();
}

void MMSDialog::preInit(){
  selectionBoxFld1->setCurrentIndex(1);
  selectionBoxFld2->setCurrentIndex(1);
  selectionBoxFld3->setCurrentIndex(1);
}

MMSDialog::~MMSDialog(){
  for (int i = 0; i < CH_NUM; i++){
    delete type1_values[i]; delete type2_values[i]; delete type3_values[i];
  }
  if (field1) delete field1; if (field2) delete field2; if (field3) delete field3;
  if (type1) delete type1; if (type2) delete type2; if (type3) delete type3;
  if (dsc_ch)      delete[] dsc_ch;
  if (peak_values) delete[] peak_values;
  if (mean_values) delete[] mean_values;
  if (sig_values)  delete[] sig_values;
  if (rms_values)  delete[] rms_values;
}

void MMSDialog::enableProceedButton()
{
  if (!runEdit->text().isEmpty() && (type1 || type2 || type3)){
    proceedButton->setEnabled(true);
    cancelButton->setDefault(false);
    proceedButton->setDefault(true);
  } else{
    proceedButton->setEnabled(false);
    proceedButton->setDefault(false);
    cancelButton->setDefault(true);
  }

  if (type1){
    selectionBoxLgc1->setEnabled(true);
    selectionBoxFld2->setEnabled(true);
    selectionBoxTyp2->setEnabled(true);
  } else {
    selectionBoxLgc1->setEnabled(false);
    selectionBoxFld2->setEnabled(false);
    selectionBoxTyp2->setEnabled(false);
  }

  if (type2){
    selectionBoxLgc2->setEnabled(true);
    selectionBoxFld3->setEnabled(true);
    selectionBoxTyp3->setEnabled(true);
  } else {
    selectionBoxLgc2->setEnabled(false);
    selectionBoxFld3->setEnabled(false);
    selectionBoxTyp3->setEnabled(false);
  }
}

void MMSDialog::setField1(int f_id){
  if (field1){ delete field1; field1 = 0; }
  if (type1) { delete type1;  type1  = 0; }
  if (f_id){
    QComboBox *box = qobject_cast<QComboBox *>(sender());
    if (box) field1 = new QString(box->itemText(f_id));
  }
  enableProceedButton();

  selectionBoxTyp1->clear();
  switch (f_id){
    case 0: return;
    case 1: for (int i = 0; i < N_types1; ++i) selectionBoxTyp1->insertItem(i, QString::fromAscii(types1[i])); break;
    case 2:
    case 3:
    case 4: for (int i = 0; i < N_types2; ++i) selectionBoxTyp1->insertItem(i, QString::fromAscii(types2[i])); break;
    case 5: for (int i = 0; i < N_types3; ++i) selectionBoxTyp1->insertItem(i, QString::fromAscii(types3[i])); break;
    default: return;
  }
}

void MMSDialog::setField2(int f_id){
  if (field2){ delete field2; field2 = 0; }
  if (type2) { delete type2;  type2  = 0; }
  if (f_id){
    QComboBox *box = qobject_cast<QComboBox *>(sender());
    if (box) field2 = new QString(box->itemText(f_id));
  }
  enableProceedButton();

  selectionBoxTyp2->clear();
  switch (f_id){
    case 0: return;
    case 1: for (int i = 0; i < N_types1; ++i) selectionBoxTyp2->insertItem(i, QString::fromAscii(types1[i])); break;
    case 2:
    case 3:
    case 4: for (int i = 0; i < N_types2; ++i) selectionBoxTyp2->insertItem(i, QString::fromAscii(types2[i])); break;
    case 5: for (int i = 0; i < N_types3; ++i) selectionBoxTyp2->insertItem(i, QString::fromAscii(types3[i])); break;
    default: return;
  }
}

void MMSDialog::setField3(int f_id){
  if (field3){ delete field3; field3 = 0; }
  if (type3) { delete type3;  type3  = 0; }
  if (f_id){
    QComboBox *box = qobject_cast<QComboBox *>(sender());
    if (box) field3 = new QString(box->itemText(f_id));
  }
  enableProceedButton();

  selectionBoxTyp3->clear();
  switch (f_id){
    case 0: return;
    case 1: for (int i = 0; i < N_types1; ++i) selectionBoxTyp3->insertItem(i, QString::fromAscii(types1[i])); break;
    case 2:
    case 3:
    case 4: for (int i = 0; i < N_types2; ++i) selectionBoxTyp3->insertItem(i, QString::fromAscii(types2[i])); break;
    case 5: for (int i = 0; i < N_types3; ++i) selectionBoxTyp3->insertItem(i, QString::fromAscii(types3[i])); break;
    default: return;
  }
}

void MMSDialog::setType1(const QString & text){
  if (type1){ delete type1; type1 = 0; }
  if (!text.isEmpty()) type1 = new QString(text);
  enableProceedButton();
}

void MMSDialog::setType2(const QString & text){
  if (type2){ delete type2; type2 = 0; }
  if (!text.isEmpty()) type2 = new QString(text);
  enableProceedButton();
}

void MMSDialog::setType3(const QString & text){
  if (type3){ delete type3; type3 = 0; }
  if (!text.isEmpty()) type3 = new QString(text);
  enableProceedButton();
}

void MMSDialog::showChannels()
{
  int t_used[3] = {1, 1, 1};
  run = runEdit->text().toInt();

  if (field1){
    if (field1->isEmpty()) { delete field1; field1 = 0; t_used[0] = 0; }
  } else t_used[0] = 0;
  if (field2){
    if (field2->isEmpty()) { delete field2; field2 = 0; t_used[1] = 0; }
  } else t_used[1] = 0;
  if (field3){
    if (field3->isEmpty()) { delete field3; field3 = 0; t_used[2] = 0; }
  } else t_used[2] = 0;
  if (type1){
    if (type1->isEmpty()) { delete type1; type1 = 0; t_used[0] = 0; }
  } else t_used[0] = 0;
  if (type2){
    if (type2->isEmpty()) { delete type2; type2 = 0; t_used[1] = 0; }
  } else t_used[1] = 0;
  if (type3){
    if (type3->isEmpty()) { delete type3; type3 = 0; t_used[2] = 0; }
  } else t_used[2] = 0;
  if ((!t_used[0] && !t_used[1] && !t_used[2]) || !run) return;

  for (int i = 0; i < CH_NUM; ++i) dsc_ch[i] = 0;

  /*
  // bad/off channels
  sprintf(the_query, "SELECT \"BadChannelsList\", \"OffChannelsList\" FROM \"LabenPrecalibDecodingQuality\" WHERE \"RunNumber\"=%d ", run);
  pg_dbi::get()->query(BX_PRECALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  int ch;
  QRegExp regExpParser(",");
  QRegExp regExpPattern("[{}]");
  QStringList bad_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(0, 0)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);
  QStringList off_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(0, 1)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);

  // disconnected channels
  sprintf(the_query, "SELECT \"ChannelID\",\"RunNumber\",\"Disconnected\" FROM \"DisconnectedPmts\" \
                      WHERE \"RunNumber\"<%d ORDER BY \"ChannelID\",\"RunNumber\" ", run+1);
  pg_dbi::get()->query(BX_CALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  char ch_state[1024];
  for (int i = 0; i < lines; ++i){
    ch = atoi(pg_dbi::get()->query_value(i, 0));
    sprintf(ch_state, pg_dbi::get()->query_value(i, 2));
    if (!strcmp(ch_state, "t")) dsc_ch[ch-1] = 1;
    else dsc_ch[ch-1] = 0;
  }

  // channels' params
  sprintf(the_query, "SELECT \"ChannelID\", \"ChargePeak\", \"ChargeMean\", \"ChargeSigma\", \"ChargeRms\" FROM \"NeutrinoPmtCalibration\" \
                      WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", run);
  pg_dbi::get()->query(BX_CALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  for (int i = 0; i < lines && i < CH_NUM; ++i){
    peak_values[i] = atof(pg_dbi::get()->query_value(i, 1));
    mean_values[i] = atof(pg_dbi::get()->query_value(i, 2));
    sig_values[i]  = atof(pg_dbi::get()->query_value(i, 3));
    rms_values[i]  = atof(pg_dbi::get()->query_value(i, 4));
  }
  */

  // channels' status
  char f1[1024], f2[1024], f3[1024];
  if (t_used[0]) sprintf(f1, field1->toLocal8Bit().constData());
  else sprintf(f1, "ChannelID");
  if (t_used[1]) sprintf(f2, field2->toLocal8Bit().constData());
  else sprintf(f2, "ChannelID");
  if (t_used[2]) sprintf(f3, field3->toLocal8Bit().constData());
  else sprintf(f3, "ChannelID");

  sprintf(the_query, "SELECT \"ChannelID\", \"%s\", \"%s\", \"%s\" FROM \"MuonChannelsProperties\" \
                      WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", f1, f2, f3, run);
  pg_dbi::get()->query(BX_CALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  QRegExp regExpPattern("[{}]");
  for (int i = 0; i < lines && i < CH_NUM; ++i){
    if (t_used[0]) sprintf(type1_values[i], (QString::fromAscii(pg_dbi::get()->query_value(i, 1)).remove(regExpPattern)).toLocal8Bit().constData());
    if (t_used[1]) sprintf(type2_values[i], (QString::fromAscii(pg_dbi::get()->query_value(i, 2)).remove(regExpPattern)).toLocal8Bit().constData());
    if (t_used[2]) sprintf(type3_values[i], (QString::fromAscii(pg_dbi::get()->query_value(i, 3)).remove(regExpPattern)).toLocal8Bit().constData());
  }

  queryContent->hide();
  queryContent->clearContents();
  queryContent->setRowCount(0);

  mainWindow::get()->dispatchMsg(str->sprintf(" >>> Muon Channels Status: run %d\n status:", run), 0);
  if (t_used[0]) mainWindow::get()->dispatchMsg(str->sprintf(" %s::%s", f1, type1->toLocal8Bit().constData()), 0);
  if (t_used[1]) mainWindow::get()->dispatchMsg(str->sprintf(" %s %s::%s", logic[selectionBoxLgc1->currentIndex()], f2, type2->toLocal8Bit().constData()), 0);
  if (t_used[2]) mainWindow::get()->dispatchMsg(str->sprintf(" %s %s::%s", logic[selectionBoxLgc2->currentIndex()], f3, type3->toLocal8Bit().constData()), 0);
  //printf("\n   lg  peak sigma  mean   rms\n");
  mainWindow::get()->dispatchMsg(str->sprintf("\n   lg\n"), 0);

  int mask, selection_found, m1, m2, m3, total_found = 0, rows_to_show = 0;
  QTableWidgetItem *tableItems0, *tableItems1, *tableItems2, *tableItems3;
  QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);

  for (int i = 0; i < lines; ++i){
    /*
    selection_found  = 0;                              // exclude disconneced, off and bad channels
    if (dsc_ch[i])  continue;
    for (int j = 0; j < bad_ch_list.size(); ++j)
      if (bad_ch_list.at(j).toInt() == i+1){ selection_found = 1; break; }
    if (selection_found) continue;
    for (int j = 0; j < off_ch_list.size(); ++j)
      if (off_ch_list.at(j).toInt() == i+1){ selection_found = 1; break; }
    if (selection_found) continue; */

    // search
    mask = selection_found = m1 = m2 = m3 = 0;
    if (t_used[0]){
      mask += 100;
      if (QString::fromAscii(type1_values[i]).contains((*type1), Qt::CaseInsensitive)) m1 = 1;
    }
    if (t_used[1]){
      mask += 100;
      mask += (selectionBoxLgc1->currentIndex() + 1)*10;
      if (QString::fromAscii(type2_values[i]).contains((*type2), Qt::CaseInsensitive)) m2 = 1;
    }
    if (t_used[2]){
      mask += 100;
      mask += (selectionBoxLgc2->currentIndex() + 1);
      if (QString::fromAscii(type3_values[i]).contains((*type3), Qt::CaseInsensitive)) m3 = 1;
    }

    //printf("mask: %d, m1: %d, m2: %d, m3: %d\n", mask, m1, m2, m3);
    switch (mask){
      case 0: break;
      case 100: if (m1) selection_found = 1; break;
      case 210: if (m1 || m2) selection_found = 1; break;
      case 220: if (m1 && m2) selection_found = 1; break;
      case 230: if (m1 && !m2) selection_found = 1; break;
      case 311: if (m1 || m2 || m3) selection_found = 1; break;
      case 312: if (m1 || m2 && m3) selection_found = 1; break;
      case 313: if (m1 || m2 && !m3) selection_found = 1; break;
      case 321: if (m1 && m2 || m3) selection_found = 1; break;
      case 322: if (m1 && m2 && m3) selection_found = 1; break;
      case 323: if (m1 && m2 && !m3) selection_found = 1; break;
      case 331: if (m1 && !m2 || m3) selection_found = 1; break;
      case 332: if (m1 && !m2 && m3) selection_found = 1; break;
      case 333: if (m1 && !m2 && !m3) selection_found = 1; break;
      default: mainWindow::get()->dispatchMsg(str->sprintf("Wrong logic management, please report this >> BUG << to developers!\n"), 0); return;
    }

    // print out
    if (selection_found){
      ++total_found;
      //printf(" %4d % 5.1f % 5.1f % 5.1f % 5.1f\n", i+1, peak_values[i], sig_values[i], mean_values[i], rms_values[i]);
      mainWindow::get()->dispatchMsg(str->sprintf(" %4d\n", i+1), 0);
    }
    else if (!checkBox->isChecked()) continue;   // leave only matching selections

    queryContent->setRowCount(rows_to_show+1);

    tableItems0 = new QTableWidgetItem(QString::number(i+1));
    tableItems0->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems0->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    if (checkBox->isChecked() && selection_found){
      tableItems0->setBackground(bg_col); tableItems0->setForeground(fg_col);
    }
    queryContent->setItem(rows_to_show, 0, tableItems0);

    if (t_used[0]){
      tableItems1 = new QTableWidgetItem(QString::fromAscii(type1_values[i]));
      tableItems1->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      tableItems1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      if (checkBox->isChecked() && selection_found){
	tableItems1->setBackground(bg_col); tableItems1->setForeground(fg_col);
      }
      queryContent->setItem(rows_to_show, 1, tableItems1);
    }

    if (t_used[1]){
      tableItems2 = new QTableWidgetItem(QString::fromAscii(type2_values[i]));
      tableItems2->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      tableItems2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      if (checkBox->isChecked() && selection_found){
	tableItems2->setBackground(bg_col); tableItems2->setForeground(fg_col);
      }
      queryContent->setItem(rows_to_show, 2, tableItems2);
    }

    if (t_used[2]){
      tableItems3 = new QTableWidgetItem(QString::fromAscii(type3_values[i]));
      tableItems3->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      tableItems3->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      if (checkBox->isChecked() && selection_found){
	tableItems3->setBackground(bg_col); tableItems3->setForeground(fg_col);
      }
      queryContent->setItem(rows_to_show, 3, tableItems3);
    }

    ++rows_to_show;
  }
  mainWindow::get()->dispatchMsg(str->sprintf(" total: %d\n", total_found), 0);

  /*
  if (prnCheckBox->isChecked()){
    total_found = 0;
    printf (" disconnected channels: ");
    for (int i = 0; i < CH_NUM; ++i){
      if (dsc_ch[i]){
	printf ("%d ", i+1);
	++total_found;
      }
    }
    mainWindow::get()->dispatchMsg(str->sprintf("\n total: %d\n", total_found);
  } */

  QStringList tableLabels; //tableLabels.clear();
  for (int i = 0; i < rows_to_show; ++i) tableLabels.append(tr("#"));
  queryContent->setVerticalHeaderLabels(tableLabels);
  tableLabels.clear();
  tableLabels.append(tr("   lg   "));
  if (t_used[0]) tableLabels.append(tr("%1::%2").arg(*field1).arg(*type1));
  else tableLabels.append(tr("                "));
  if (t_used[1]) tableLabels.append(tr("%1::%2").arg(*field2).arg(*type2));
  else tableLabels.append(tr("                "));
  if (t_used[2]) tableLabels.append(tr("%1::%2").arg(*field3).arg(*type3));
  else tableLabels.append(tr("                "));
  queryContent->setHorizontalHeaderLabels(tableLabels);

  queryContent->resizeColumnsToContents();
  //queryContent->adjustSize();
  queryContent->update();
  queryContent->show();

  update();
  fflush(stdout);
}

