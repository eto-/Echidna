#ifndef _BX_DB_MON_MCH_MON_DIALOG_H
#define _BX_DB_MON_MCH_MON_DIALOG_H

#include "dbMon.h"
#include <QDialog>

//QT_BEGIN_NAMESPACE
class QLabel;
class QLineEdit;
class QButtonGroup;
class QRadioButton;
class QPushButton;
class QTableWidget;
class QCheckBox;
//QT_END_NAMESPACE


class MCMDialog : public QDialog
{
  Q_OBJECT

  public:
    MCMDialog(QWidget *parent = 0);
    ~MCMDialog();

  public slots:
    void enableProceedButton();
    void showChannels();

  private:
    int step, run1, run2, grpn;
    int lines, columns;
    char the_query[4096];
    int *run_n; 	//, *off_ch, *sum_ch, *bad_ls, *off_ls;
    float *tt1_rt, *tt2_rt, *mean_rt;
    QString *str;

    QTableWidget *queryContent;

    QLabel *run1Label, *run2Label, *stepLabel, *meanLabel;
    QLineEdit *run1Edit, *run2Edit, *stepEdit, *meanEdit;

    QButtonGroup *methRadGrp, *typeRadGrp;
    QRadioButton *prevRadio, *averRadio, *tt1Radio, *tt2Radio;
    QCheckBox    *meanRadio;
    QPushButton *proceedButton, *cancelButton;
};

#endif
