#include <QtGui>
#include <math.h>

#include "mainWindow.h"
#include "dmpmondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QTableWidget>
#include <QStringList>
#include <QCheckBox>

const int   DMPDialog::HV_SET_MIN        = 100;
const int   DMPDialog::DHV_MIN           = 10;
const float DMPDialog::DARK_NOISE_THRES  = 2.0; // kHz
const float DMPDialog::BAD_SHAPE_THRES   = 10.;
const float DMPDialog::LOW_CHARGE_THRES  = 20.;
const float DMPDialog::HIGH_CHARGE_THRES = 30.;
const char  DMPDialog::tbl_clm[][1024]   = {"PRECALIB", "DEAD_IN_NU", "LOW_GAIN", "HIGH_GAIN", "HOT_IN_NU", "BAD_SHAPE", "LOW_CHARGE", "HIGH_CHARGE", "BAD_HV"};

// *********************** outline window

DMPDialog::DMPDialog(QWidget *parent) : QDialog(parent), tagged(0)
{
    str = new QString();

    tagged = new int*[CH_NUM];
    for (int i = 0; i < CH_NUM; ++i) tagged[i] = new int[N_clm];
    for (int i = 0; i < CH_NUM; ++i){
	multStatus[i] = new char[4096]; chStatus[i] = new char[4096]; chBaseS[i] = new char[4096]; chPeakS[i] = new char[4096]; tmStatus[i] = new char[4096];
    }
    dsc_ch = new int[CH_NUM];

    crate  = new int[CH_NUM];
    hv_set = new int[CH_NUM];
    hv_mon = new int[CH_NUM];
    feb    = new float[CH_NUM];
    lbn    = new float[CH_NUM];
    hvb    = new float[CH_NUM];

    peak_values = new float[CH_NUM];
    mean_values = new float[CH_NUM];
    sig_values  = new float[CH_NUM];
    rms_values  = new float[CH_NUM];
    toff_values = new float[CH_NUM];
    tsig_values = new float[CH_NUM];
    dkns_values = new float[CH_NUM];

    QRegExp regExp("[0-9]{0,6}");
    runLabel = new QLabel(tr("&RUN:"), this);
    runEdit  = new QLineEdit(this);
    runLabel->setBuddy(runEdit);
    //runEdit->setText(tr("7733"));
    runEdit->setText(tr("%1").arg(mainWindow::current_run()));
    runEdit->setValidator(new QRegExpValidator(regExp, this));

    connect(runEdit, SIGNAL(textChanged(const QString &)), this, SLOT(enableProceedButton()));

    proceedButton = new QPushButton(tr("&Find Channels"), this);
    proceedButton->setEnabled(true);
    dumpButton = new QPushButton(tr("&Dump to DB"), this);
    dumpButton->setEnabled(false);
    dumpButton->setDefault(false);
    cancelButton = new QPushButton(tr("&Cancel"), this);
    cancelButton->setDefault(true);
    cancelButton->setEnabled(true);

    connect(proceedButton, SIGNAL(clicked()), this, SLOT(showChannels()));
    connect(dumpButton, SIGNAL(clicked()), this, SLOT(dumpChannels()));
    connect(cancelButton, SIGNAL(clicked()), this, SLOT(close()));

    queryContent = new QTableWidget(0, N_clm+1, this);
    queryContent->hide();

    QHBoxLayout *controlLayout = new QHBoxLayout;
    controlLayout->addWidget(runLabel);
    controlLayout->addWidget(runEdit);
    controlLayout->addWidget(proceedButton);
    controlLayout->addWidget(dumpButton);
    controlLayout->insertStretch(5, 1);
    controlLayout->addWidget(cancelButton);

    QVBoxLayout *tableLayout = new QVBoxLayout;
    tableLayout->addWidget(queryContent, 1);

    QVBoxLayout *mainLayout = new QVBoxLayout;
    mainLayout->addLayout(controlLayout);
    mainLayout->addLayout(tableLayout, 1);
    setLayout(mainLayout);

    setMinimumWidth(controlLayout->sizeHint().width());
    //setFixedWidth(sizeHint().width());

    setWindowTitle(tr("Laben Broken Channels"));

    mapChannels();
}

DMPDialog::~DMPDialog(){
    if (tagged){
	for (int i = 0; i < CH_NUM; ++i){
	    delete[] tagged[i];
	}
	delete[] tagged;
    }
    for (int i = 0; i < CH_NUM; i++){
	delete multStatus[i]; delete chStatus[i]; delete chBaseS[i]; delete chPeakS[i]; delete tmStatus[i];
    }
    if (dsc_ch)      delete[] dsc_ch;
    if (crate)       delete[] crate;
    if (feb)         delete[] feb;
    if (lbn)         delete[] lbn;
    if (hvb)         delete[] hvb;
    if (peak_values) delete[] peak_values;
    if (mean_values) delete[] mean_values;
    if (sig_values)  delete[] sig_values;
    if (rms_values)  delete[] rms_values;
    if (toff_values) delete[] toff_values;
    if (tsig_values) delete[] tsig_values;
    if (dkns_values) delete[] dkns_values;
    if (hv_set)      delete[] hv_set;
    if (hv_mon)      delete[] hv_mon;
}

void DMPDialog::enableProceedButton()
{
    if (!runEdit->text().isEmpty()){
	proceedButton->setEnabled(true);
	cancelButton->setDefault(false);
	proceedButton->setDefault(true);
    } else{
	proceedButton->setEnabled(false);
	proceedButton->setDefault(false);
	cancelButton->setDefault(true);
    }
    dumpButton->setEnabled(false);
}

void DMPDialog::mapChannels()
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
    pg_dbi::get()->query(1, the_query);
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


void DMPDialog::dumpChannels()
{
    int dump2db;
    char reason[2048];

    sprintf(the_query, "SELECT \"RunNumber\" FROM \"BrokenChannels\" WHERE \"RunNumber\"=%d ", run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    if (lines){
	int ret = QMessageBox::question(this, tr("RUN %1 DB DUMP").arg(run), tr("DB entries for run %1 already exist. "
		                                 "Proceed with dump (existing entries will be overwritten)?").arg(run),
		                                 QMessageBox::Ok | QMessageBox::Cancel, QMessageBox::Cancel);
	switch (ret) {
	    case QMessageBox::Ok:
		mainWindow::get()->dispatchMsg(str->sprintf(" Deleting old DB entries for run %d...", run), 0); fflush(stdout);
		break;
	    case QMessageBox::Cancel:
		dumpButton->setEnabled(false);
		return;
	    default:
		break;
	}

	sprintf(the_query, "DELETE FROM \"BrokenChannels\" WHERE \"RunNumber\"=%d ", run);
	pg_dbi::get()->query(BX_CALIB, the_query, 1);
    }

    //for (int i = 0; i < CH_NUM; ++i){
    //    if (i < 780 || i > 790) continue;
    //    mainWindow::get()->dispatchMsg(str->sprintf("ch %4d crate %2d feb % 4.2f laben % 4.2f hvb % 4.2f\n", i+1, crate[i], feb[i], lbn[i], hvb[i]),0);
    //}

    sprintf(reason, "./laben_broken_channels_%06d.txt", run);
    FILE *fascii_out = 0;
    if (!(fascii_out = fopen(reason, "w"))) {
	mainWindow::get()->dispatchMsg(str->sprintf(" Can't open %s, dump to ASCII file aborted.\n", reason), 0); fflush(stdout);
    }

    for (int i = 0; i < CH_NUM; ++i){
	reason[0] = '\0';
	dump2db = 0;
	for (int j = 0; j < N_clm; ++j){
	    if (tagged[i][j]){
		if ((j == 1) && (tagged[i][0])) continue; // write out PRECALIB && DEAD_IN_NU as PRECALIB only
		if (dump2db) strcat(reason, ", ");
		strcat(reason, tbl_clm[j]);
		dump2db = 1;
	    }
	}
	if (!dump2db) continue;

	sprintf(the_query, "INSERT INTO \"BrokenChannels\" (\"RunNumber\", \"Reason\", \"ChannelID\", \"Crate\", \"Feb\", \"Lbnb\", \"Hvb\", \
                 	    \"ChargeBaseStatus\", \"ChargePeakStatus\", \"ChargeStatus\", \"TimingStatus\", \"LabenMultiplicity\", \
			    \"ChargePeak\", \"ChargeSigma\", \"ChargeMean\", \"ChargeRms\", \"TimeOffset\", \"TimeSigma\", \"HvSet\", \"HvMon\", \"DarkRate\") \
			    VALUES (%d, \'%s\', %d, %d, %5.2f, %5.2f, %5.2f, \
			    \'%s\', \'%s\', \'%s\', \'%s\', \'%s\', \
			    %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %6.3f)",
		            run, reason, i+1, crate[i], feb[i], lbn[i], hvb[i],
		            chBaseS[i], chPeakS[i], chStatus[i], tmStatus[i], multStatus[i],
		            peak_values[i], sig_values[i], mean_values[i], rms_values[i], toff_values[i], tsig_values[i], hv_set[i], hv_mon[i], dkns_values[i]);
	pg_dbi::get()->query(BX_CALIB, the_query, 1);

        reason[0] = '\0';
        if (tagged[i][0]) strcat(reason, "1"); else strcat(reason, "2");
	if (fascii_out) fprintf(fascii_out, "%d %s\n", i+1, reason);
    }

    if (fascii_out) fclose(fascii_out);

    dumpButton->setEnabled(false);
    mainWindow::get()->dispatchMsg(str->sprintf(" dump completed.\n"), 0); fflush(stdout);
}


void DMPDialog::showChannels()
{
    int ch, i, j, lv_ch;
    float HV_AVG, TM_AVG;

    run = runEdit->text().toInt();
    if (!run) return;

    for (i = 0; i < CH_NUM; ++i){
	dsc_ch[i] = 0;
	for (j = 0; j < N_clm; ++j) tagged[i][j] = 0;
    }

    // disconnected channels
    //char ch_state[BUFSIZE];
    sprintf(the_query, "SELECT \"ChannelID\", \"RunNumber\", \"Disconnected\" FROM \"DisconnectedPmts\" \
	                WHERE \"RunNumber\"<%d ORDER BY \"ChannelID\", \"RunNumber\" ", run+1);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

//     mainWindow::get()->dispatchMsg(str->sprintf("  *** mapping before\n"), 0);
//     for (ch = 0; ch < CH_NUM; ++ch){
// 	if (ch < 786 || ch > 790) continue;
// 	mainWindow::get()->dispatchMsg(str->sprintf("ch %4d crate %2d feb % 4.2f laben % 4.2f hvb % 4.2f\n", ch+1, crate[ch], feb[ch], lbn[ch], hvb[ch]),0);
//     }

    for (i = 0; i < lines; ++i){
	ch = atoi(pg_dbi::get()->query_value(i, 0));
	if (!strncmp(pg_dbi::get()->query_value(i, 2), "t", 1)) {
	    dsc_ch[ch-1] = 1;
	    //mainWindow::get()->dispatchMsg(str->sprintf("  %d diconnected\n", ch-1), 0);
	} else dsc_ch[ch-1] = 0;  // obligatory!!!
    }

//     mainWindow::get()->dispatchMsg(str->sprintf("  *** mapping after\n"), 0);
//     for (ch = 0; ch < CH_NUM; ++ch){
// 	if (ch < 786 || ch > 790) continue;
// 	mainWindow::get()->dispatchMsg(str->sprintf("ch %4d crate %2d feb % 4.2f laben % 4.2f hvb % 4.2f\n", ch+1, crate[ch], feb[ch], lbn[ch], hvb[ch]),0);
//     }

    mapChannels();  // patching memory...

    // bad/off channels
    QRegExp regExpParser(",");
    QRegExp regExpPattern("[{}]");
    sprintf(the_query, "SELECT \"BadChannelsList\", \"OffChannelsList\" FROM \"LabenPrecalibDecodingQuality\" WHERE \"RunNumber\"=%d ", run);
    pg_dbi::get()->query(BX_PRECALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    QStringList bad_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(0, 0)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);
    QStringList off_ch_list = (QString::fromAscii(pg_dbi::get()->query_value(0, 1)).remove(regExpPattern)).split(regExpParser, QString::SkipEmptyParts);

    // search for the closest calibration run
    int calib_run = run+1; lines = 0;
    sprintf(the_query, "SELECT \"RunNumber\" FROM \"LaserPmtCalibration\" ORDER BY \"RunNumber\" ");
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    while (lines && calib_run > run)
	calib_run = QString::fromAscii(pg_dbi::get()->query_value(--lines, 0)).toInt();
    //printf("calib_run: %d, lines: %d\n", calib_run, lines);
    if (calib_run < 5003) calib_run = -1;

    // search for the closest HV info run
    int  found = 0, hv_run = run > first_hv_run ? run+1 : -1;
    char hv_dmp_reason[1024] = "";
    while (hv_run > first_hv_run && !found){
	--hv_run;
	sprintf(the_query, "SELECT \"Reason\" FROM \"HighVoltages\" WHERE \"Run\"=%d ", hv_run);
	pg_dbi::get()->query(BX_SLOW, the_query);
	pg_dbi::get()->query_dimensions(lines, columns);

	for (i = 0; i < lines; ++i){
	    sprintf(hv_dmp_reason, pg_dbi::get()->query_value(i, 0));
	    if (QString::fromAscii(hv_dmp_reason).contains(("RGDMP"), Qt::CaseInsensitive)){ found = 1;  break; }
	}
        //printf("hv_run: %d, lines: %d\n", hv_run, lines); fflush(stdout);
    }
    if (hv_run == first_hv_run) hv_run = -1;

    // channels' statuses
    sprintf(the_query, "SELECT \"ChannelID\", \"Multiplicity\", \"ChargeStatus\", \"ChargeBaseStatus\", \"ChargePeakStatus\", \"TimingStatus\" \
	                FROM \"LabenChannelsProperties\" WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ",  run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (i = 0; i < lines && i < CH_NUM; ++i){
	sprintf(multStatus[i], (QString::fromAscii(pg_dbi::get()->query_value(i, 1)).remove(regExpPattern)).toLocal8Bit().constData());
	sprintf(chStatus[i],   (QString::fromAscii(pg_dbi::get()->query_value(i, 2)).remove(regExpPattern)).toLocal8Bit().constData());
	sprintf(chBaseS[i],    (QString::fromAscii(pg_dbi::get()->query_value(i, 3)).remove(regExpPattern)).toLocal8Bit().constData());
	sprintf(chPeakS[i],    (QString::fromAscii(pg_dbi::get()->query_value(i, 4)).remove(regExpPattern)).toLocal8Bit().constData());
	sprintf(tmStatus[i],   (QString::fromAscii(pg_dbi::get()->query_value(i, 5)).remove(regExpPattern)).toLocal8Bit().constData());
    }

    // channels' params
    sprintf(the_query, "SELECT \"ChannelID\", \"ChargePeak\", \"ChargeMean\", \"ChargeSigma\", \"ChargeRms\" FROM \"NeutrinoPmtCalibration\" \
	                WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (i = 0; i < lines && i < CH_NUM; ++i){
	peak_values[i] = atof(pg_dbi::get()->query_value(i, 1));
	mean_values[i] = atof(pg_dbi::get()->query_value(i, 2));
	sig_values[i]  = atof(pg_dbi::get()->query_value(i, 3));
	rms_values[i]  = atof(pg_dbi::get()->query_value(i, 4));
    }

    // time values
    sprintf(the_query, "SELECT \"ChannelID\", \"TimeOffset\", \"TimeSigma\" FROM \"LaserPmtCalibration\" \
	                WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", calib_run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    TM_AVG = 0.; lv_ch = 0;
    for (i = 0; i < lines && i < CH_NUM; ++i){
	toff_values[i] = atof(pg_dbi::get()->query_value(i, 1));
	tsig_values[i] = atof(pg_dbi::get()->query_value(i, 2));
	TM_AVG += toff_values[i]; ++lv_ch;
    }
    if (lv_ch) TM_AVG /= static_cast<float>(lv_ch);

    // PMT's dark noise values
    sprintf(the_query, "SELECT \"ChannelID\", \"DarkNoise\" FROM \"InnerPmtsDarkRate\" \
	                WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", calib_run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    for (i = 0; i < lines && i < CH_NUM; ++i){
	dkns_values[i] = atof(pg_dbi::get()->query_value(i, 1));
    }

    // HV info
    if (hv_run > 0){
	sprintf(the_query, "SELECT \"Crate\", \"Board\", \"Reason\", \
	              	   \"V0SET1\",  \"VMON1\",  \"V0SET2\",  \"VMON2\",  \"V0SET3\",  \"VMON3\",  \"V0SET4\",  \"VMON4\", \
			   \"V0SET5\",  \"VMON5\",  \"V0SET6\",  \"VMON6\",  \"V0SET7\",  \"VMON7\",  \"V0SET8\",  \"VMON8\", \
			   \"V0SET9\",  \"VMON9\",  \"V0SET10\", \"VMON10\", \"V0SET11\", \"VMON11\", \"V0SET12\", \"VMON12\", \
			   \"V0SET13\", \"VMON13\", \"V0SET14\", \"VMON14\", \"V0SET15\", \"VMON15\", \"V0SET16\", \"VMON16\", \
			   \"V0SET17\", \"VMON17\", \"V0SET18\", \"VMON18\", \"V0SET19\", \"VMON19\", \"V0SET20\", \"VMON20\", \
			   \"V0SET21\", \"VMON21\", \"V0SET22\", \"VMON22\", \"V0SET23\", \"VMON23\", \"V0SET24\", \"VMON24\"  \
			   FROM \"HighVoltages\" WHERE \"Run\"=%d ORDER BY \"Time\" DESC", hv_run);
	pg_dbi::get()->query(BX_SLOW, the_query);
	pg_dbi::get()->query_dimensions(lines, columns);
	mainWindow::get()->dispatchMsg(str->sprintf(" >>> Laben Broken Channels: run %d, calibration run %d, HV info run %d\n", run, calib_run, hv_run), 0);
    } else {
	mainWindow::get()->dispatchMsg(str->sprintf(" >>> Laben Broken Channels: run %d, calibration run %d, HV info run NOT FOUND in DB.\n", run, calib_run), 0);
    }

    HV_AVG = 0.; lv_ch = 0;
    int crt, brd;
    for (i = 0; i < CH_NUM; ++i){
	hv_set[i] = hv_mon[i] = -1;
	if (dsc_ch[i]) continue;

	if (hv_run > 0){
	    crt = crate[i];
	    brd = static_cast<int>(hvb[i]);
	    ch  = static_cast<int>(101.*hvb[i] - 101.*brd);

	    for (j = 0; j < lines; ++j){
		if( (crt == atoi(pg_dbi::get()->query_value(j, 0))) && (brd == atoi(pg_dbi::get()->query_value(j, 1))) &&
				   (QString::fromAscii(pg_dbi::get()->query_value(j, 2)).contains(("RGDMP"), Qt::CaseInsensitive)) ){
		    hv_set[i] = atoi(pg_dbi::get()->query_value(j, 2*ch+1));
		    hv_mon[i] = atoi(pg_dbi::get()->query_value(j, 2*ch+2));
		    HV_AVG += hv_set[i]; ++lv_ch;
		    break;
		}
	    }
	}
    }
    if (lv_ch) HV_AVG /= static_cast<float>(lv_ch);


    queryContent->hide();
    queryContent->clearContents();
    queryContent->setRowCount(0);

    mainWindow::get()->dispatchMsg(str->sprintf("   lg   peak  sigma   mean    rms tm_off tm_sig drk_ns hv_set hv_mon   reason(s)\n"), 0);

    int selection_found, total_found = 0, rows_to_show = 0;
    QTableWidgetItem *tableItems0, *tableItems1, *tableItems2, *tableItems3, *tableItems4, *tableItems5, *tableItems6, *tableItems7, *tableItems8, *tableItems9;
    QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);

    for (i = 0; i < CH_NUM; ++i){
	selection_found = 0;
	if (dsc_ch[i]) continue;    // exclude disconneced

        // 1 problem in precalibration
	for (j = 0; j < bad_ch_list.size(); ++j) if (bad_ch_list.at(j).toInt() == i+1)                  { selection_found = 1; tagged[i][0] = 1; }
	for (j = 0; j < off_ch_list.size(); ++j) if (off_ch_list.at(j).toInt() == i+1)                  { selection_found = 1; tagged[i][0] = 1; }

        // 2 dead in neutrino
	if (QString::fromAscii(multStatus[i]).contains(("dead_in_neutrino"), Qt::CaseInsensitive))          { selection_found = 1; tagged[i][1] = 1; }

        // 3 low gain
	if (QString::fromAscii(chStatus[i]).contains(("low_gain"),  Qt::CaseInsensitive))                   { selection_found = 1; tagged[i][2] = 1; }

        // 4 high gain
	if (QString::fromAscii(chStatus[i]).contains(("high_gain"), Qt::CaseInsensitive))                   { selection_found = 1; tagged[i][3] = 1; }

        // 5 hot in neutrino
	if (dkns_values[i] > 0. && dkns_values[i] < DARK_NOISE_THRES &&
            QString::fromAscii(multStatus[i]).contains(("hot_in_neutrino"), Qt::CaseInsensitive) &&
            !QString::fromAscii(multStatus[i]).contains(("retriggering_in_neutrino"), Qt::CaseInsensitive)) { selection_found = 1; tagged[i][4] = 1; }

        // 6 bad shape of charge peak
	if (/*peak_values[i] > 0. && */fabs(peak_values[i] - mean_values[i]) > BAD_SHAPE_THRES)             { selection_found = 1; tagged[i][5] = 1; }

        // 7 low charge
	if (/*peak_values[i] > 0. && */peak_values[i] < LOW_CHARGE_THRES)                                   { selection_found = 1; tagged[i][6] = 1; }

        // 8 high charge
	if (peak_values[i] > HIGH_CHARGE_THRES)                                                             { selection_found = 1; tagged[i][7] = 1; }

        // 9 high charge
	if (hv_set[i]*hv_mon[i] > 1 && hv_set[i] > HV_SET_MIN && abs(hv_set[i]-hv_mon[i]) > DHV_MIN)        { selection_found = 1; tagged[i][8] = 1; }

        //printf("ch %4d crate %2d hvb % 4.2f hv_set % .0f hv_mon % .0f\n", i+1, crate[i], hvb[i], hv_set[i], hv_mon[i]); fflush(stdout);
	if (selection_found){
	    ++total_found;

            // print out
	    mainWindow::get()->dispatchMsg(str->sprintf(" %4d % 6.1f % 6.1f % 6.1f % 6.1f % 6.1f % 6.1f % 6.2f % 6d % 6d :", i+1, peak_values[i], sig_values[i], mean_values[i], rms_values[i], toff_values[i], tsig_values[i], dkns_values[i], hv_set[i], hv_mon[i]), 0);
	    for (j = 0; j < N_clm; ++j){
		if (j == 1 && tagged[i][0]) continue;
		if (tagged[i][j]) mainWindow::get()->dispatchMsg(str->sprintf(" %s :", tbl_clm[j]), 0); fflush(stdout);
	    }
	    mainWindow::get()->dispatchMsg(str->sprintf("\n"), 0);
	}
	else continue;   // leave only matching selections

	queryContent->setRowCount(rows_to_show+1);

	tableItems0 = new QTableWidgetItem(QString::number(i+1));
	tableItems0->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
	tableItems0->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 0, tableItems0);

	tableItems1 = new QTableWidgetItem();
	if (tagged[i][0]){
            //tableItems1->setBackground(bg_col); tableItems1->setForeground(fg_col);
	    tableItems1->setText("X");
	} else {
	    tableItems1->setText("");
	}
	tableItems1->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems1->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 1, tableItems1);

	tableItems2 = new QTableWidgetItem();
	if (tagged[i][1]){
            //tableItems2->setBackground(bg_col); tableItems2->setForeground(fg_col);
	    tableItems2->setText("X");
	} else {
	    tableItems2->setText("");
	}
	tableItems2->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems2->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 2, tableItems2);

	tableItems3 = new QTableWidgetItem();
	if (tagged[i][2]){
            //tableItems3->setBackground(bg_col); tableItems3->setForeground(fg_col);
	    tableItems3->setText("X");
	} else {
	    tableItems3->setText("");
	}
	tableItems3->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems3->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 3, tableItems3);

	tableItems4 = new QTableWidgetItem();
	if (tagged[i][3]){
            //tableItems4->setBackground(bg_col); tableItems4->setForeground(fg_col);
	    tableItems4->setText("X");
	} else {
	    tableItems4->setText("");
	}
	tableItems4->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems4->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 4, tableItems4);

	tableItems5 = new QTableWidgetItem();
	if (tagged[i][4]){
            //tableItems5->setBackground(bg_col); tableItems5->setForeground(fg_col);
	    tableItems5->setText("X");
	} else {
	    tableItems5->setText("");
	}
	tableItems5->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems5->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 5, tableItems5);

	tableItems6 = new QTableWidgetItem();
	if (tagged[i][5]){
            //tableItems6->setBackground(bg_col); tableItems6->setForeground(fg_col);
	    tableItems6->setText("X");
	} else {
	    tableItems6->setText("");
	}
	tableItems6->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems6->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 6, tableItems6);

	tableItems7 = new QTableWidgetItem();
	if (tagged[i][6]){
            //tableItems7->setBackground(bg_col); tableItems7->setForeground(fg_col);
	    tableItems7->setText("X");
	} else {
	    tableItems7->setText("");
	}
	tableItems7->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems7->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 7, tableItems7);

	tableItems8 = new QTableWidgetItem();
	if (tagged[i][7]){
            //tableItems8->setBackground(bg_col); tableItems8->setForeground(fg_col);
	    tableItems8->setText("X");
	} else {
	    tableItems8->setText("");
	}
	tableItems8->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems8->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 8, tableItems8);

	tableItems9 = new QTableWidgetItem();
	if (tagged[i][8]){
            //tableItems9->setBackground(bg_col); tableItems9->setForeground(fg_col);
	    tableItems9->setText("X");
	} else {
	    tableItems9->setText("");
	}
	tableItems9->setTextAlignment(Qt::AlignHCenter | Qt::AlignVCenter);
	tableItems9->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
	queryContent->setItem(rows_to_show, 9, tableItems9);

	++rows_to_show;
    }
    mainWindow::get()->dispatchMsg(str->sprintf(" total: %d\n", total_found), 0);

    QStringList tableLabels; //tableLabels.clear();
    for (i = 0; i < rows_to_show; ++i) tableLabels.append(tr("#"));
    queryContent->setVerticalHeaderLabels(tableLabels);
    tableLabels.clear();
    tableLabels.append(tr("   lg   "));
    for (i = 0; i < N_clm; ++i) tableLabels.append(tr("%1").arg(tbl_clm[i]));
    queryContent->setHorizontalHeaderLabels(tableLabels);

    queryContent->resizeColumnsToContents();
    //queryContent->adjustSize();
    queryContent->update();
    queryContent->show();

    dumpButton->setEnabled(true);

    update();
    fflush(stdout);
}

/*    sprintf(the_query, "SELECT \"VMON1\",  \"V0SET1\",  \"VMON2\",  \"V0SET2\",  \"VMON3\",  \"V0SET3\",  \"VMON4\",  \"V0SET4\",  \"VMON5\",  \"V0SET5\",  \"VMON6\",  \"V0SET6\",  \
	         	       \"VMON7\",  \"V0SET7\",  \"VMON8\",  \"V0SET8\",  \"VMON9\",  \"V0SET9\",  \"VMON10\", \"V0SET10\", \"VMON11\", \"V0SET11\", \"VMON12\", \"V0SET12\", \
			       \"VMON13\", \"V0SET13\", \"VMON14\", \"V0SET14\", \"VMON15\", \"V0SET15\", \"VMON16\", \"V0SET16\", \"VMON17\", \"V0SET17\", \"VMON18\", \"V0SET18\", \
			       \"VMON19\", \"V0SET19\", \"VMON20\", \"V0SET20\", \"VMON21\", \"V0SET21\", \"VMON22\", \"V0SET22\", \"VMON23\", \"V0SET23\", \"VMON24\", \"V0SET24\"  \
		        	FROM \"HighVoltages\" WHERE \"Run\"=%d AND \"Crate\"=%d AND \"Board\"=%d ORDER BY \"Run\" ", hv_run, crate[i], brd);*/
