#include <QtGui>
#include <math.h>

#include "mainWindow.h"
#include "hvmondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QTableWidget>
#include <QStringList>
#include <QCheckBox>

const char  HVDialog::tbl_clm[][1024]   = {"RACK", "BOARD", "CHANNEL", "HV_SET", "HV_MON", "DISCONNECTED"};

// *********************** outline window

HVDialog::HVDialog(QWidget *parent) : QDialog(parent), tagged(0)
{
  str = new QString();

  tagged = new int*[CH_NUM];
  for (int i = 0; i < CH_NUM; ++i) tagged[i] = new int[N_clm];
  dsc_ch = new int[CH_NUM];

  crate  = new int[CH_NUM];
  feb    = new float[CH_NUM];
  lbn    = new float[CH_NUM];
  hvb    = new float[CH_NUM];
  hv_set = new int[CH_NUM];
  hv_mon = new int[CH_NUM];
  hv_dbl = new int[CH_NUM];

  QRegExp regExp("[0-9]{0,6}");
  runLabel = new QLabel(tr("&RUN:"), this);
  runEdit  = new QLineEdit(this);
  runLabel->setBuddy(runEdit);
  runEdit->setText(tr("%1").arg(mainWindow::current_run()));
  runEdit->setValidator(new QRegExpValidator(regExp, this));

  connect(runEdit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));

  QRegExp regExp1("[0-9]{0,4}");
  dhvLabel = new QLabel(tr("  Voltage &diff >"), this);
  sr1Label = new QLabel(tr("V "), this);
  dhvEdit  = new QLineEdit(this);
  dhvLabel->setBuddy(dhvEdit);
  dhvEdit->setText(tr("10"));
  dhvEdit->setValidator(new QRegExpValidator(regExp1, this));
  dhvEdit->setFixedWidth((int)(0.4*dhvEdit->sizeHint().width()));

  hvsLabel = new QLabel(tr("  HV &set >"), this);
  hvsEdit  = new QLineEdit(this);
  sr2Label = new QLabel(tr("V   "), this);
  hvsLabel->setBuddy(hvsEdit);
  hvsEdit->setText(tr("100"));
  hvsEdit->setValidator(new QRegExpValidator(regExp1, this));
  hvsEdit->setFixedWidth((int)(0.4*hvsEdit->sizeHint().width()));

  proceedButton = new QPushButton(tr("&Find Channels"), this);
  proceedButton->setEnabled(true);

  connect(proceedButton, SIGNAL(clicked()), this, SLOT(showChannels()));

  queryContent = new QTableWidget(0, N_clm+1, this);
  queryContent->hide();

  prnLabel = new QLabel(tr("print out disconnected"), this);
  prnCheckBox = new QCheckBox(this);
  prnCheckBox->setChecked(true);

  QHBoxLayout *controlLayout = new QHBoxLayout;
  controlLayout->addWidget(runLabel);
  controlLayout->addWidget(runEdit);
  controlLayout->addWidget(dhvLabel);
  controlLayout->addWidget(dhvEdit);
  controlLayout->addWidget(sr1Label);
  controlLayout->addWidget(hvsLabel);
  controlLayout->addWidget(hvsEdit);
  controlLayout->addWidget(sr2Label);
  controlLayout->insertStretch(8, 1);
  controlLayout->addWidget(proceedButton);

  QVBoxLayout *tableLayout = new QVBoxLayout;
  tableLayout->addWidget(queryContent, 1);

  QHBoxLayout *boxLayout = new QHBoxLayout;
  //boxLayout->addWidget(checkBox);
  //boxLayout->addWidget(boxLabel);
  boxLayout->insertStretch(1, 1);
  boxLayout->addWidget(prnLabel);
  boxLayout->addWidget(prnCheckBox);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addLayout(controlLayout);
  mainLayout->addLayout(tableLayout, 1);
  mainLayout->addLayout(boxLayout);
  setLayout(mainLayout);

  setMinimumWidth(controlLayout->sizeHint().width());
  //setFixedWidth(sizeHint().width());

  setWindowTitle(tr("HV Channels Monitor"));

  mapChannels();
}

HVDialog::~HVDialog(){
  if (tagged){
    for (int i = 0; i < CH_NUM; ++i){
      delete[] tagged[i];
    }
    delete[] tagged;
  }
  if (dsc_ch)      delete[] dsc_ch;
  if (crate)       delete[] crate;
  if (feb)         delete[] feb;
  if (lbn)         delete[] lbn;
  if (hvb)         delete[] hvb;
  if (hv_set)      delete[] hv_set;
  if (hv_mon)      delete[] hv_mon;
  if (hv_dbl)      delete[] hv_dbl;
}

void HVDialog::enableProceedButton()
{
  if (!runEdit->text().isEmpty()){
    proceedButton->setEnabled(true);
    proceedButton->setDefault(true);
  } else{
    proceedButton->setEnabled(false);
    proceedButton->setDefault(false);
  }
}

void HVDialog::mapChannels()
{
  int crt, lbnb, lbnc;
  int lbnbd[CH_NUM], lbnch[CH_NUM];

  for (int ch = 0; ch < CH_NUM; ++ch){
    crt  = static_cast<int>(floor((ch)/160.) + 1.);            // crate
    lbnb = static_cast<int>(floor((ch-(crt-1)*160.)/8.) + 1);  // Laben board
    lbnc = (ch+1-(crt-1)*160-(lbnb-1)*8);                      // Laben channel

    lbnbd[ch] = lbnb; lbnch[ch] = lbnc;

    crate[ch] = crt;
    lbn[ch] = lbnb + static_cast<float>(lbnc)/100.;
    feb[ch] = hvb[ch] = 0.;
    //printf("ch %4d laben %4.2f \n", ch+1, lbn[ch]);
  }

  sprintf(the_query, "SELECT \"RackID\", \"FEBoard\", \"FEChannel\", \"LabenBoard\", \"LabenChannel\", \"HVBoard\", \"HVChannel\"  \
                      FROM  \"CableMapping\" WHERE \"ProfileID\"=1 ORDER BY \"RackID\" ");
  pg_dbi::get()->query(BX_GEOMETRY, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  for (int i = 0; i < lines; ++i){
    crt  = atoi(pg_dbi::get()->query_value(i, 0));
    lbnb = atoi(pg_dbi::get()->query_value(i, 3));
    lbnc = atoi(pg_dbi::get()->query_value(i, 4));
    for (int ch = 0; ch < CH_NUM; ++ch){
      if (crate[ch] != crt || lbnbd[ch] != lbnb || lbnch[ch] != lbnc) continue;
      //crate[ch] = atoi(pg_dbi::get()->query_value(i, 0));
      feb[ch]   = atof(pg_dbi::get()->query_value(i, 1)) + atof(pg_dbi::get()->query_value(i, 2))/100.;
      hvb[ch]   = atof(pg_dbi::get()->query_value(i, 5)) + atof(pg_dbi::get()->query_value(i, 6))/100.;
      //printf("ch %4d crate %2d feb % 4.2f laben % 4.2f hvb % 4.2f\n", ch+1, crate[ch], feb[ch], lbn[ch], hvb[ch]); fflush(stdout);
      break;
    }
  }
}

void HVDialog::showChannels()
{
  run = runEdit->text().toInt();
  if (!run) return;

  dhv = dhvEdit->text().toInt();
  hvs = hvsEdit->text().toInt();
  //printf("Voltage: %d, min hv_set: %d\n", dhv, hvs); fflush(stdout);
  for (int i = 0; i < CH_NUM; ++i){
    dsc_ch[i] = 0;
    for (int j = 0; j < N_clm; ++j) tagged[i][j] = 0;
  }

  // disconnected channels
  sprintf(the_query, "SELECT \"ChannelID\",\"RunNumber\",\"Disconnected\" FROM \"DisconnectedPmts\" \
                      WHERE \"RunNumber\"<%d ORDER BY \"ChannelID\",\"RunNumber\" ", run+1);
  pg_dbi::get()->query(BX_CALIB, the_query);
  pg_dbi::get()->query_dimensions(lines, columns);

  int ch;
  char ch_state[1024];
  for (int i = 0; i < lines; ++i){
    ch = atoi(pg_dbi::get()->query_value(i, 0));
    sprintf(ch_state, pg_dbi::get()->query_value(i, 2));
    if (!strcmp(ch_state, "t")) dsc_ch[ch-1] = 1;
    else dsc_ch[ch-1] = 0;
  }

  // search for the closest HV info run
  int  found = 0, hv_run = run > first_hv_run ? run+1 : -1;
  char hv_dmp_reason[1024] = "";
  while (hv_run > first_hv_run && !found){
    --hv_run;
    sprintf(the_query, "SELECT \"Reason\" FROM \"HighVoltages\" WHERE \"Run\"=%d ", hv_run);
    pg_dbi::get()->query(BX_SLOW, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (int i = 0; i < lines; i++){
      sprintf(hv_dmp_reason, pg_dbi::get()->query_value(i, 0));
      if (QString::fromAscii(hv_dmp_reason).contains(("RGDMP"), Qt::CaseInsensitive)){ found = 1;  break; }
    }
    //printf("hv_run: %d, lines: %d\n", hv_run, lines); fflush(stdout);
  }
  if (hv_run == first_hv_run) hv_run = -1;

  if (hv_run > 0){
    sprintf(the_query, "SELECT \"Crate\", \"Board\", \
	                \"V0SET1\",  \"VMON1\",  \"V0SET2\",  \"VMON2\",  \"V0SET3\",  \"VMON3\",  \"V0SET4\",  \"VMON4\", \
			\"V0SET5\",  \"VMON5\",  \"V0SET6\",  \"VMON6\",  \"V0SET7\",  \"VMON7\",  \"V0SET8\",  \"VMON8\", \
			\"V0SET9\",  \"VMON9\",  \"V0SET10\", \"VMON10\", \"V0SET11\", \"VMON11\", \"V0SET12\", \"VMON12\", \
			\"V0SET13\", \"VMON13\", \"V0SET14\", \"VMON14\", \"V0SET15\", \"VMON15\", \"V0SET16\", \"VMON16\", \
			\"V0SET17\", \"VMON17\", \"V0SET18\", \"VMON18\", \"V0SET19\", \"VMON19\", \"V0SET20\", \"VMON20\", \
			\"V0SET21\", \"VMON21\", \"V0SET22\", \"VMON22\", \"V0SET23\", \"VMON23\", \"V0SET24\", \"VMON24\", \"Reason\"  \
			FROM \"HighVoltages\" WHERE \"Run\"=%d ORDER BY \"Time\" DESC", hv_run);
    pg_dbi::get()->query(BX_SLOW, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    mainWindow::get()->dispatchMsg(str->sprintf(" >>> HV Channels Status: run %d, HV info run %d\n", run, hv_run), 0);
  } else {
    mainWindow::get()->dispatchMsg(str->sprintf(" >>> HV Channels Status: run %d, HV info run NOT FOUND in DB. Skipping further processing.\n", run), 0);
    return;
  }

  queryContent->hide();
  queryContent->clearContents();
  queryContent->setRowCount(0);

  mainWindow::get()->dispatchMsg(str->sprintf(" threshold values: abs(hv_set-hv_mon) > %d V, hv_set > %d V\n", dhv, hvs), 0);
  mainWindow::get()->dispatchMsg(str->sprintf("   lg  crate  board  hv_ch hv_set hv_mon discon\n"), 0);

  int crt, brd, selection_found, total_found = 0, rows_to_show = 0;
  QTableWidgetItem *tableItems0, *tableItems1, *tableItems2, *tableItems3, *tableItems4, *tableItems5, *tableItems6;
  QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);

  for (int i = 0; i < CH_NUM; ++i){
    selection_found = 0;

    if (!prnCheckBox->isChecked() && dsc_ch[i]) continue;

    // HV settings, if exist
    crt = crate[i];
    brd = static_cast<int>(hvb[i]);
    ch  = static_cast<int>(101.*hvb[i] - 101.*brd);

    hv_set[i] = hv_mon[i] = -1; hv_dbl[i] = 0;
    for (int j = 0; j < lines; ++j){
      if( (crt == atoi(pg_dbi::get()->query_value(j, 0))) && (brd == atoi(pg_dbi::get()->query_value(j, 1))) &&
	  (QString::fromAscii(pg_dbi::get()->query_value(j, 50)).contains(("RGDMP"), Qt::CaseInsensitive)) ){
	if (hv_set[i]*hv_mon[i] == 1){
	  hv_set[i] = atoi(pg_dbi::get()->query_value(j, 2*ch));
	  hv_mon[i] = atoi(pg_dbi::get()->query_value(j, 2*ch+1));
	  //break;
	} else {
	  hv_dbl[i] = 1;
	  break;
	}
      }
      else continue;
    }
    if (hv_set[i] < hvs) continue;
    if (hv_set[i]*hv_mon[i] != 1 && abs(hv_set[i]-hv_mon[i]) >= dhv) selection_found = 1;

    if (selection_found){
      ++total_found;
      mainWindow::get()->dispatchMsg(str->sprintf(" %4d   %4d   %4d   %4d   %4d   %4d   %2d", i+1, crt, brd, ch, hv_set[i], hv_mon[i], dsc_ch[i]), 0);
      if (hv_dbl[i]) mainWindow::get()->dispatchMsg(str->sprintf("     WARNING: entry's not unique, showing the latest"), 0);
      mainWindow::get()->dispatchMsg(str->sprintf("\n"), 0);
    } else continue;

    queryContent->setRowCount(rows_to_show+1);

    tableItems0 = new QTableWidgetItem(QString::number(i+1));
    tableItems0->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
    tableItems0->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 0, tableItems0);

    tableItems1 = new QTableWidgetItem(QString::number(crt));
    tableItems1->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 1, tableItems1);

    tableItems2 = new QTableWidgetItem(QString::number(brd));
    tableItems2->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 2, tableItems2);

    tableItems3 = new QTableWidgetItem(QString::number(ch));
    tableItems3->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems3->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 3, tableItems3);

    tableItems4 = new QTableWidgetItem(QString::number(hv_set[i]));
    tableItems4->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems4->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 4, tableItems4);

    tableItems5 = new QTableWidgetItem(QString::number(hv_mon[i]));
    tableItems5->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems5->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 5, tableItems5);

    tableItems6 = new QTableWidgetItem();
    if (dsc_ch[i]){
      //tableItems6->setBackground(bg_col); tableItems6->setForeground(fg_col);
      tableItems6->setText("X");
    } else {
      tableItems6->setText("");
    }
    tableItems6->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
    tableItems6->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
    queryContent->setItem(rows_to_show, 6, tableItems6);

    ++rows_to_show;
  }
  mainWindow::get()->dispatchMsg(str->sprintf(" total: %d\n", total_found), 0);

  QStringList tableLabels; //tableLabels.clear();
  for (int i = 0; i < rows_to_show; ++i) tableLabels.append(tr("#"));
  queryContent->setVerticalHeaderLabels(tableLabels);
  tableLabels.clear();
  tableLabels.append(tr("   lg   "));
  for (int i = 0; i < N_clm; ++i) tableLabels.append(tr("%1").arg(tbl_clm[i]));
  queryContent->setHorizontalHeaderLabels(tableLabels);

  queryContent->resizeColumnsToContents();
  //  queryContent->adjustSize();
  queryContent->update();
  queryContent->show();

  update();
  fflush(stdout);
}
