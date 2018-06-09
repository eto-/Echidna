#ifndef _BX_DB_MON_DMPMU_MON_DIALOG_H
#define _BX_DB_MON_DMPMU_MON_DIALOG_H

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



class DMPMuDialog : public QDialog
{
  Q_OBJECT

  public:
    DMPMuDialog(QWidget *parent = 0);
    ~DMPMuDialog();

  public slots:
    void enableProceedButton();
    void showChannels();
    void dumpChannels();

  private:
    void mapChannels();

    static const int   N_clm        = 7;
    static const int   first_hv_run = 16910;
    static const int   PROF_ID;
    static const int   DHV_MIN;
    static const int   HV_SET_MIN;
    static const float CHARGE_MIN_THRES;
    static const float CHARGE_MAX_THRES;
    static const float CHARGE_SIG_THRES;
    static const float TIME_THRES;
    static const float TIME_SIG_THRES;
    static const char  tbl_clm[][1024];

    int lines, columns;
    char the_query[4096];
    QString *str;

    int run;
    int **tagged;
    int *dsc_ch, *crate, *ppl, *hid, *hv_set, *hv_mon;
    float *qtc, *tdc, *hvb;
    float *peak_values, *sig_values, *toff_values, *tsig_values, *dkns_values;
    char  *multStatus[CH_NUM_MU];

    QLineEdit *runEdit;
    QLabel *runLabel, *boxLabel, *prnLabel; // *selectionLabel1, *selectionLabel2, *selectionLabel3;
    QCheckBox *checkBox, *prnCheckBox;

    QTableWidget *queryContent;
    QPushButton *proceedButton, *dumpButton, *cancelButton;
};

#endif
