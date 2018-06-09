#ifndef _BX_DB_MON_DB_DIALOG_H
#define _BX_DB_MON_DB_DIALOG_H

#include "dbMon.h"
#include <QDialog>

//QT_BEGIN_NAMESPACE
class QTabWidget;
class QListWidget;
class QListWidgetItem;
class QTableWidget;
class QTableWidgetItem;
class QDialogButtonBox;
class QLabel;
class QCheckBox;
//QT_END_NAMESPACE

class infoTab : public QWidget
{
  Q_OBJECT

  public:
    infoTab(const int, QWidget *parent = 0);
    ~infoTab(){ destroy_arrays(); }

    void setDB(const int);

  public slots:
    void showTableFields(QListWidgetItem*);
    void prnt_state(int state) { prnt = state; }

  private:
    void destroy_arrays();

    int db_id, n_tables, prnt;
    int lines, columns;
    char the_query[4096];
    QString *str;

    int  *tbl_n_fields;
    char ***tbl_fields_name;
    char ***tbl_fields_type;

    QListWidget  *tablesListBox;
    //QListWidgetItem  *listItems;
    QTableWidget *tableContent;
    //QTableWidgetItem *tableItemsName, *tableItemsType;
};

class queryTab : public QWidget
{
  Q_OBJECT

  public:
    queryTab(const int, QWidget *parent = 0);
    ~queryTab(){ }

    void setDB(const int id) { db_id = id; };

  public slots:
    //void showTableFields(QListWidgetItem*);
    //void prnt_state(int state) { prnt = state; }

  private:
    //void destroy_arrays();

    int db_id;
    int lines, columns;
    char the_query[4096];
    QString *str;

/*    int  *tbl_n_fields;
    char ***tbl_fields_name;
    char ***tbl_fields_type;
    QListWidget  *tablesListBox;*/
    //QListWidgetItem  *listItems;
    //QTableWidget *tableContent;
    //QTableWidgetItem *tableItemsName, *tableItemsType;
};

class DBDialog : public QDialog
{
  Q_OBJECT

  public:
    DBDialog(const int, QWidget *parent = 0);

  public slots:
    void setDBId();

  private:
    int db_id;
    QString *str;

    QTabWidget *tabWidget;
    QLabel *prnLabel;
    QCheckBox *prnCheckBox;
    QDialogButtonBox *buttonBox;

    infoTab  *tab1;
    queryTab *tab2;
};

#endif
