#include <QtGui>
#include <math.h>

#include "mainWindow.h"
#include "dmpmumondialog.h"

#include <QLabel>
#include <QLineEdit>
#include <QComboBox>
#include <QPushButton>
#include <QTableWidget>
#include <QStringList>
#include <QCheckBox>

const int   DMPMuDialog::PROF_ID          = 19;
const int   DMPMuDialog::DHV_MIN          = 10;
const int   DMPMuDialog::HV_SET_MIN       = 100;
const float DMPMuDialog::CHARGE_MIN_THRES = 5.;
const float DMPMuDialog::CHARGE_MAX_THRES = 30.;
const float DMPMuDialog::CHARGE_SIG_THRES = 13.;
const float DMPMuDialog::TIME_THRES       = 10.;
const float DMPMuDialog::TIME_SIG_THRES   = 3.;
const char  DMPMuDialog::tbl_clm[][1024]   = {"DEAD_IN_NU", "LOW_EFF_NU", "LOW_CHARGE", "HIGH_CHARGE", "BAD_SHAPE", "BAD_TIMING", "BAD_HV"};

// *********************** outline window

DMPMuDialog::DMPMuDialog(QWidget *parent) : QDialog(parent), tagged(0)
{
    str = new QString();

    tagged = new int*[CH_NUM_MU];
    for (int i = 0; i < CH_NUM_MU; ++i) tagged[i] = new int[N_clm];
    for (int i = 0; i < CH_NUM_MU; ++i){
	multStatus[i] = new char[4096];
    }
    dsc_ch = new int[CH_NUM_MU];

    crate  = new int[CH_NUM_MU];
    ppl    = new int[CH_NUM_MU];
    hid    = new int[CH_NUM_MU];
    hv_set = new int[CH_NUM_MU];
    hv_mon = new int[CH_NUM_MU];
    qtc    = new float[CH_NUM_MU];
    tdc    = new float[CH_NUM_MU];
    hvb    = new float[CH_NUM_MU];

    peak_values = new float[CH_NUM_MU];
    sig_values  = new float[CH_NUM_MU];
    toff_values = new float[CH_NUM_MU];
    tsig_values = new float[CH_NUM_MU];
    dkns_values = new float[CH_NUM_MU];

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

    setWindowTitle(tr("Muon Broken Channels"));

    mapChannels();
}

DMPMuDialog::~DMPMuDialog(){
    if (tagged){
	for (int i = 0; i < CH_NUM_MU; ++i) delete[] tagged[i];
	delete[] tagged;
    }
    for (int i = 0; i < CH_NUM_MU; i++){
	delete multStatus[i];
    }
    if (dsc_ch)      delete[] dsc_ch;
    if (crate)       delete[] crate;
    if (ppl)         delete[] ppl;
    if (hid)         delete[] hid;
    if (qtc)         delete[] qtc;
    if (tdc)         delete[] tdc;
    if (hvb)         delete[] hvb;
    if (peak_values) delete[] peak_values;
    if (sig_values)  delete[] sig_values;
    if (toff_values) delete[] toff_values;
    if (tsig_values) delete[] tsig_values;
    if (dkns_values) delete[] dkns_values;
    if (hv_set)      delete[] hv_set;
    if (hv_mon)      delete[] hv_mon;
}

void DMPMuDialog::enableProceedButton()
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

void DMPMuDialog::mapChannels()
{
    for (int ch = 0; ch < CH_NUM_MU; ++ch){
	crate[ch] = ppl[ch] = hid[ch] = 0;
	qtc[ch] = tdc[ch] = hvb[ch] = 0.;
	if ((ch < 112) && (ch%16))               { /* ppl[ch] = ch - int(ch/16);      */ crate[ch] = 15; }
	if ((ch > 128) && (ch < 238) && (ch%16)) { /* ppl[ch] = ch - int(ch/16) - 15; */ crate[ch] = 15; }
    }

    sprintf(the_query, "SELECT \"ChannelID\", \"HoleLabel\" FROM  \"MuonHolesMapping\" WHERE \"ProfileID\"=%d ORDER BY \"ChannelID\" ", PROF_ID);
    pg_dbi::get()->query(BX_GEOMETRY, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (int i = 0; i < lines; ++i) hid[atoi(pg_dbi::get()->query_value(i, 0)) - 3001]  = atoi(pg_dbi::get()->query_value(i, 1));

    sprintf(the_query, "SELECT \"HoleLabel\", \"PatchPanelLabel\", \"QTCBoard\", \"QTCChannel\", \
	                \"TDCBoard\", \"TDCChip\", \"TDCChannel\", \"HVBoard\", \"HVChannel\" \
			FROM  \"MuonCableMapping\" WHERE \"ProfileID\"=%d ORDER BY \"PatchPanelLabel\" ", PROF_ID);
    pg_dbi::get()->query(BX_GEOMETRY, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    int hidv;
    for (int i = 0; i < lines; ++i){
	hidv  = atoi(pg_dbi::get()->query_value(i, 0));
	for (int ch = 0; ch < CH_NUM_MU; ++ch){
	    if (!crate[ch] || hid[ch] != hidv) continue;
	    ppl[ch]   = atoi(pg_dbi::get()->query_value(i, 1));
	    qtc[ch]   = atof(pg_dbi::get()->query_value(i, 2)) + atof(pg_dbi::get()->query_value(i, 3))/100.;
	    tdc[ch]   = atof(pg_dbi::get()->query_value(i, 4)) + (32.*atof(pg_dbi::get()->query_value(i, 5)) + atof(pg_dbi::get()->query_value(i, 6)))/100.;
	    hvb[ch]   = atof(pg_dbi::get()->query_value(i, 7)) + atof(pg_dbi::get()->query_value(i, 8))/100.;
	    //printf("ch %4d cr %2d ppl %3d hid % 3d qtc % 5.2f tdc % 5.2f hvb % 5.2f\n", ch+3001, crate[ch], ppl[ch], hid[ch], qtc[ch], tdc[ch], hvb[ch]); fflush(stdout);
	    break;
	}
    }
}

void DMPMuDialog::dumpChannels()
{
    int dump2db;
    char reason[2048];

    sprintf(the_query, "SELECT \"RunNumber\" FROM \"BrokenMuonChannels\" WHERE \"RunNumber\"=%d ", run);
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

	sprintf(the_query, "DELETE FROM \"BrokenMuonChannels\" WHERE \"RunNumber\"=%d ", run);
	pg_dbi::get()->query(BX_CALIB, the_query, 1);
    }

    sprintf(reason, "./muon_broken_channels_%06d.txt", run);
    FILE *fascii_out = 0;
    if (!(fascii_out = fopen(reason, "w"))) {
	mainWindow::get()->dispatchMsg(str->sprintf(" Can't open %s, dump to ASCII file aborted.\n", reason), 0); fflush(stdout);
    }

    for (int i = 0; i < CH_NUM_MU; ++i){
	reason[0] = '\0';
	dump2db = 0;
	for (int j = 0; j < N_clm; ++j){
	    if (tagged[i][j]){
		//if ((j == 1) && (tagged[i][0])) continue; // write out PRECALIB && DEAD_IN_NU as PRECALIB only
		if (dump2db) strcat(reason, ", ");
		strcat(reason, tbl_clm[j]);
		dump2db = 1;
	    }
	}
	if (!dump2db) continue;

	sprintf(the_query, "INSERT INTO \"BrokenMuonChannels\" (\"ChannelID\", \"RunNumber\", \"Reason\", \"Mch\", \"Crate\", \"PPL\", \"QTC\", \"TDC\", \"Hvb\", \
		            \"Multiplicity\", \"ChargePeak\", \"ChargeSigma\", \"TimeOffset\", \"TimeSigma\", \"HvSet\", \"HvMon\", \"DarkRate\") \
			    VALUES (%d, %d, \'%s\', %d, %d, %d, %5.2f, %5.2f, %5.2f, \
			    \'%s\', %0.2f, %0.2f, %0.2f, %0.2f, %d, %d, %6.3f)",
		            3001 + i, run, reason, i, crate[i], ppl[i], qtc[i], tdc[i], hvb[i],
		            multStatus[i], peak_values[i], sig_values[i], toff_values[i], tsig_values[i], hv_set[i], hv_mon[i], dkns_values[i]);
	pg_dbi::get()->query(BX_CALIB, the_query, 1);

        //reason[0] = '\0';
        //if (tagged[i][0]) strcat(reason, "1"); else strcat(reason, "2");
	if (fascii_out) fprintf(fascii_out, "%d %s\n", i+3001, reason);
    }

    if (fascii_out) fclose(fascii_out);

    dumpButton->setEnabled(false);
    mainWindow::get()->dispatchMsg(str->sprintf(" dump completed.\n"), 0);
}

void DMPMuDialog::showChannels()
{
    int ch, lv_ch;
    float HV_AVG, TM_AVG;

    run = runEdit->text().toInt();
    if (!run) return;

    for (int i = 0; i < CH_NUM_MU; ++i){
	dsc_ch[i] = 0;
	for (int j = 0; j < N_clm; ++j) tagged[i][j] = 0;
    }

    // disconnected channels
    sprintf(the_query, "SELECT \"ChannelID\", \"RunNumber\", \"Disconnected\" FROM \"DisconnectedPmts\" \
	                WHERE \"RunNumber\"<%d ORDER BY \"ChannelID\", \"RunNumber\" ", run+1);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (int i = 0; i < lines; ++i){
	ch = atoi(pg_dbi::get()->query_value(i, 0)) - 3001;
	if (ch < 0 || ch > CH_NUM_MU) continue;
	if (!strcmp(pg_dbi::get()->query_value(i, 2), "t")) dsc_ch[ch] = 1;
	else dsc_ch[ch] = 0; // obligatory!!!
    }

    // empty channles
    for (int i = 0; i < CH_NUM_MU; ++i) if (!ppl[i]) dsc_ch[i] = 1;

    mapChannels();  // memory patch...

    // search for the closest calibration run
    int calib_run = run+1; lines = 0;
    sprintf(the_query, "SELECT \"RunNumber\" FROM \"MuonPmtCalibration\" ORDER BY \"RunNumber\" ");
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    while (lines && calib_run > run)
	calib_run = QString::fromAscii(pg_dbi::get()->query_value(--lines, 0)).toInt();
    //printf("calib_run: %d, lines: %d\n", calib_run, lines);
    if (calib_run < 5003) calib_run = -1;

    // search for the closest HV run
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

    // channels' statuses
    sprintf(the_query, "SELECT \"ChannelID\", \"Multiplicity\" \
	                FROM \"MuonChannelsProperties\" WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ",  run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    QRegExp regExpPattern("[{}]");
    for (int i = 0; i < lines; ++i){
	ch = atoi(pg_dbi::get()->query_value(i, 0)) - 3001;
	if (ch < 0 || ch > CH_NUM_MU) continue;
	sprintf(multStatus[ch], (QString::fromAscii(pg_dbi::get()->query_value(i, 1)).remove(regExpPattern)).toLocal8Bit().constData());
    }

    // channels' params
    sprintf(the_query, "SELECT \"ChannelID\", \"ChargePeak\", \"ChargeSigma\", \"TimeOffset\", \"TimeSigma\" FROM \"MuonPmtCalibration\" \
	                WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", calib_run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    TM_AVG = 0.; lv_ch = 0;
    for (int i = 0; i < lines; ++i){
	ch = atoi(pg_dbi::get()->query_value(i, 0)) - 3001;
	if (ch < 0 || ch > CH_NUM_MU) continue;
	peak_values[ch] = atof(pg_dbi::get()->query_value(i, 1));
	sig_values[ch]  = atof(pg_dbi::get()->query_value(i, 2));
	toff_values[ch] = atof(pg_dbi::get()->query_value(i, 3));
	tsig_values[ch] = atof(pg_dbi::get()->query_value(i, 4));
	TM_AVG += toff_values[ch]; ++lv_ch;
    }
    if (lv_ch) TM_AVG /= static_cast<float>(lv_ch);

    // PMT's dark noise values
    sprintf(the_query, "SELECT \"ChannelID\", \"DarkNoise\" FROM \"OuterPmtsDarkRate\" \
	                WHERE \"RunNumber\"=%d ORDER BY \"ChannelID\" ", calib_run);
    pg_dbi::get()->query(BX_CALIB, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);

    for (int i = 0; i < lines; ++i){
	ch = atoi(pg_dbi::get()->query_value(i, 0)) - 3001;
	if (ch < 0 || ch > CH_NUM_MU) continue;
	dkns_values[ch] = atof(pg_dbi::get()->query_value(i, 1));
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
	mainWindow::get()->dispatchMsg(str->sprintf(" >>> Muon Broken Channels: run %d, calibration run %d, HV info run %d\n", run, calib_run, hv_run), 0);
    } else {
	mainWindow::get()->dispatchMsg(str->sprintf(" >>> Muon Broken Channels: run %d, calibration run %d, HV info run NOT FOUND in DB.\n", run, calib_run), 0);
    }

    HV_AVG = 0.; lv_ch = 0;
    int crt, brd;
    for (int i = 0; i < CH_NUM_MU; ++i){
	hv_set[i] = hv_mon[i] = -1;
	if (dsc_ch[i]) continue;

	if (hv_run > 0){
	    crt = crate[i];
	    brd = static_cast<int>(hvb[i]);
	    ch  = static_cast<int>(101.*hvb[i] - 101.*brd);

	    for (int j = 0; j < lines; ++j){
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

    //mainWindow::get()->dispatchMsg(str->sprintf("HV_AVG %f, TM_AVG %f\n", HV_AVG, TM_AVG), 0);

    queryContent->hide();
    queryContent->clearContents();
    queryContent->setRowCount(0);

    mainWindow::get()->dispatchMsg(str->sprintf("   lg   peak  sigma tm_off tm_sig drk_ns hv_set hv_mon   reason(s)\n"), 0);

    int selection_found, total_found = 0, rows_to_show = 0;
    QTableWidgetItem *tableItems0, *tableItems1, *tableItems2, *tableItems3, *tableItems4, *tableItems5, *tableItems6, *tableItems7; // *tableItems8, *tableItems9;
    QBrush fg_col(Qt::white); QBrush bg_col(Qt::red);

    for (int i = 0; i < CH_NUM_MU; ++i){
	selection_found = 0;
	//printf("ppl: %d, mch: %d, dsc: %d, ch: %d\n", ppl[i], i, dsc_ch[i], 3001 + i); fflush(stdout);
	if (dsc_ch[i]) continue;    // exclude disconneced

        // 1 dead in neutrino
	if (QString::fromAscii(multStatus[i]).contains(("dead_in_neutrino"), Qt::CaseInsensitive))          { selection_found = 1; tagged[i][0] = 1; }

	// 2 low efficiency in neutrino
	if (QString::fromAscii(multStatus[i]).contains(("low_eff_in_neutrino"), Qt::CaseInsensitive))       { selection_found = 1; tagged[i][1] = 1; }

        // 3 low charge
	if (peak_values[i] < CHARGE_MIN_THRES)                                                              { selection_found = 1; tagged[i][2] = 1; }

        // 4 high charge
	if (peak_values[i] > CHARGE_MAX_THRES)                                                              { selection_found = 1; tagged[i][3] = 1; }

	// 5 bad shape
	if (sig_values[i] > CHARGE_SIG_THRES)                                                               { selection_found = 1; tagged[i][4] = 1; }

        // 6 bad timing
	if (toff_values[i] < -TIME_THRES || toff_values[i] > TIME_THRES || tsig_values[i] > TIME_SIG_THRES) { selection_found = 1; tagged[i][5] = 1; }

	// 7 bad HV
	if (hv_set[i]*hv_mon[i] > 1 && hv_set[i] > HV_SET_MIN && abs(hv_set[i]-hv_mon[i]) > DHV_MIN)        { selection_found = 1; tagged[i][6] = 1; }

	if (selection_found){
	    ++total_found;

            // print out
	    mainWindow::get()->dispatchMsg(str->sprintf(" %4d % 6.1f % 6.1f % 6.1f % 6.1f %6.2f % 6d % 6d :", i+3001, peak_values[i], sig_values[i], toff_values[i], tsig_values[i], dkns_values[i], hv_set[i], hv_mon[i]), 0);
	    for (int j = 0; j < N_clm; ++j){
		if (j == 1 && tagged[i][0]) continue;
		if (tagged[i][j]) mainWindow::get()->dispatchMsg(str->sprintf(" %s :", tbl_clm[j]), 0);
	    }
	    mainWindow::get()->dispatchMsg(str->sprintf("\n"), 0);
	}
	else continue;   // leave only matching selections

	queryContent->setRowCount(rows_to_show+1);

	tableItems0 = new QTableWidgetItem(QString::number(i+3001));
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

    dumpButton->setEnabled(true);

    update();
    fflush(stdout);
}
