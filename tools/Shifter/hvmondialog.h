#ifndef _BX_DB_MON_HV_MON_DIALOG_H
#define _BX_DB_MON_HV_MON_DIALOG_H

#include "dbMon.h"
#include <QDialog>

//QT_BEGIN_NAMESPACE
class QLabel;
class QLineEdit;
class QComboBox;
class QPushButton;
class QTableWidget;
class QString;
class QCheckBox;
//QT_END_NAMESPACE

class HVDialog : public QDialog
{
  Q_OBJECT

  public:
    HVDialog(QWidget *parent = 0);
    ~HVDialog();

  public slots:
    void enableProceedButton();
    void showChannels();
    //void dumpChannels();

  private:
    void mapChannels();

    static const int N_clm        = 6;
    static const int first_hv_run = 16910;
    static const char tbl_clm[][1024];

    int lines, columns;
    char the_query[4096];
    QString *str;

    int run, dhv, hvs;
    int **tagged;
    int *dsc_ch, *crate, *hv_set, *hv_mon, *hv_dbl;
    float *feb, *lbn, *hvb;
    //float *peak_values, *mean_values, *sig_values, *rms_values, *toff_values, *tsig_values, *dkns_values;
    //char  *multStatus[CH_NUM], *chStatus[CH_NUM], *chBaseS[CH_NUM], *chPeakS[CH_NUM], *tmStatus[CH_NUM];

    QLineEdit *runEdit, *dhvEdit, *hvsEdit;
    QLabel *runLabel, *dhvLabel, *sr1Label, *hvsLabel, *sr2Label, /*, *boxLabel, */ *prnLabel;
    QCheckBox /* *checkBox, */ *prnCheckBox;

    QTableWidget *queryContent;
    QPushButton *proceedButton; //, *dumpButton, *cancelButton;
};

#endif
