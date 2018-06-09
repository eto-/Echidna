#ifndef _BX_DB_MON_LMS_MON_DIALOG_H
#define _BX_DB_MON_LMS_MON_DIALOG_H

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

class LMSDialog : public QDialog
{
  Q_OBJECT

  public:
    LMSDialog(QWidget *parent = 0);
    ~LMSDialog();

  public slots:
    void enableProceedButton();
    void showChannels();
    void setField1(int);
    void setField2(int);
    void setField3(int);
    void setType1(const QString &);
    void setType2(const QString &);
    void setType3(const QString &);

  private:
    static const int N_conds  = 3;
    static const int N_fields = 6;  // remember add  +1 for unselected option
    static const int N_types1 = 21; // remember add  +1 for unselected option
    static const int N_types2 = 10; // remember add  +1 for unselected option
    static const int N_types3 = 9;  // remember add  +1 for unselected option

    static const char logic[][1024];
    static const char fields[][1024];
    static const char types1[][1024];
    static const char types2[][1024];
    static const char types3[][1024];

    int lines, columns;
    char the_query[4096];
    QString *str;

    int run;
    int *dsc_ch;
    float *peak_values, *mean_values, *sig_values, *rms_values;
    char *type1_values[CH_NUM], *type2_values[CH_NUM], *type3_values[CH_NUM];
    QString *field1, *field2, *field3;
    QString *type1, *type2, *type3;

    QLineEdit *runEdit;
    QLabel *runLabel, *boxLabel, *prnLabel; // *selectionLabel1, *selectionLabel2, *selectionLabel3;
    QComboBox *selectionBoxFld1, *selectionBoxFld2, *selectionBoxFld3;
    QComboBox *selectionBoxTyp1, *selectionBoxTyp2, *selectionBoxTyp3;
    QComboBox *selectionBoxLgc1, *selectionBoxLgc2;
    QCheckBox *checkBox, *prnCheckBox;

    QTableWidget *queryContent;
    QPushButton *proceedButton, *cancelButton;
};

#endif
