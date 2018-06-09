#ifndef _BX_DB_MON_DMP_MON_DIALOG_H
#define _BX_DB_MON_DMP_MON_DIALOG_H

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

class DMPDialog : public QDialog
{
  Q_OBJECT

  public:
    DMPDialog(QWidget *parent = 0);
    ~DMPDialog();

  public slots:
    void enableProceedButton();
    void showChannels();
    void dumpChannels();

  private:
    void mapChannels();

    static const int   N_clm        = 9;
    static const int   first_hv_run = 16910;
    static const int   HV_SET_MIN;
    static const int   DHV_MIN;
    static const float DARK_NOISE_THRES;
    static const float BAD_SHAPE_THRES;
    static const float LOW_CHARGE_THRES;
    static const float HIGH_CHARGE_THRES;
    static const char  tbl_clm[][1024];

    int lines, columns, length;
    char the_query[4096];
    QString *str;

    int run;
    int **tagged;
    int *dsc_ch, *crate, *hv_set, *hv_mon;
    float *feb, *lbn, *hvb;
    float *peak_values, *mean_values, *sig_values, *rms_values, *toff_values, *tsig_values, *dkns_values;
    char  *multStatus[CH_NUM], *chStatus[CH_NUM], *chBaseS[CH_NUM], *chPeakS[CH_NUM], *tmStatus[CH_NUM];

    QLineEdit *runEdit;
    QLabel *runLabel, *boxLabel, *prnLabel; // *selectionLabel1, *selectionLabel2, *selectionLabel3;
    QCheckBox *checkBox, *prnCheckBox;

    QTableWidget *queryContent;
    QPushButton *proceedButton, *dumpButton, *cancelButton;
};

#endif
