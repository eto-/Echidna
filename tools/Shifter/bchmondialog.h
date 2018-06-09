#ifndef _BX_DB_MON_BCH_MON_DIALOG_H
#define _BX_DB_MON_BCH_MON_DIALOG_H

#include "dbMon.h"
#include <QDialog>

//QT_BEGIN_NAMESPACE
class QLabel;
class QLineEdit;
class QRadioButton;
class QPushButton;
class QTableWidget;
//QT_END_NAMESPACE


class BCMDialog : public QDialog
{
  Q_OBJECT

  public:
    BCMDialog(QWidget *parent = 0);
    ~BCMDialog();

  public slots:
    void enableProceedButton();
    void showChannels();

  private:
    int step, run1, run2, pcnt;
    int lines, columns;
    char the_query[4096];
    int *run_n, *bad_ch, *off_ch, *sum_ch, *bad_ls, *off_ls;
    QString *str;

    QTableWidget *queryContent;

    QLabel *run1Label, *run2Label, *stepLabel, *pcntLabel;
    QLineEdit *run1Edit, *run2Edit, *stepEdit, *pcntEdit;

    QRadioButton *prevRadio, *averRadio;
    QPushButton *proceedButton, *cancelButton;
};

#endif
