#include <QtGui>
// #include <QTabWidget>
// #include <QLabel>
// #include <QCheckBox>
// #include <QDialogButtonBox>

#include "dbdialog.h"

// *********************** outline window

DBDialog::DBDialog(const int id, QWidget *parent)
  : QDialog(parent), db_id(id)
{
  str = new QString();

  tab1 = new infoTab(db_id, this),
  tab2 = new queryTab(db_id, this);
  //tab3 = new PermissionsTab(db_id, this);

  tabWidget = new QTabWidget;
  tabWidget->addTab(tab1, tr("DB content"));
  //tabWidget->addTab(tab2, tr("Submit query"));
  //tabWidget->addTab(tab3, tr("Permissions"));

  prnLabel = new QLabel(tr("print out"), this);
  prnCheckBox = new QCheckBox(this);
  prnCheckBox->setChecked(false);

  connect(prnCheckBox, SIGNAL(stateChanged (int)), tab1, SLOT(prnt_state (int)));

  buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok);
  buttonBox->button(QDialogButtonBox::Ok)->setAutoDefault(false);
  buttonBox->button(QDialogButtonBox::Ok)->setDefault(false);
  connect(buttonBox, SIGNAL(accepted()), this, SLOT(accept()));
  //connect(buttonBox, SIGNAL(rejected()), this, SLOT(reject()));

  QHBoxLayout *floorLayout = new QHBoxLayout;
  floorLayout->addWidget(prnCheckBox);
  floorLayout->addWidget(prnLabel);
  floorLayout->insertStretch(3, 1);
  floorLayout->addWidget(buttonBox);

  QVBoxLayout *mainLayout = new QVBoxLayout;
  mainLayout->addWidget(tabWidget, 1);
  mainLayout->addLayout(floorLayout);
  setLayout(mainLayout);

  setMinimumWidth((int)(sizeHint().width()));

  QWidget::setTabOrder(tab1, tab2);
  QWidget::setTabOrder(tab2, buttonBox);
  QWidget::setTabOrder(buttonBox, tab1);
  setWindowTitle(tr("%1").arg(db_names[db_id]));
}

void DBDialog::setDBId(){
  QAction *action = qobject_cast<QAction *>(sender());
  if (action) db_id = action->data().toInt();

  tab1->setDB(db_id);
  tab1->update();

  tab2->setDB(db_id);
  tab2->update();

  buttonBox->button(QDialogButtonBox::Ok)->setDefault(false);
  setWindowTitle(tr("%1").arg(db_names[db_id]));
}

//int DBDialog::printout() {return prnCheckBox->isChecked(); }

// *********************** tab1

infoTab::infoTab(const int id, QWidget *parent)
  : QWidget(parent), prnt(0), tbl_n_fields(0)
{
  str = new QString();

  tablesListBox = new QListWidget(this);
  connect(tablesListBox, SIGNAL(itemActivated(QListWidgetItem*)), this, SLOT(showTableFields(QListWidgetItem*)));
  connect(tablesListBox, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(showTableFields(QListWidgetItem*)));
  //connect(tablesListBox, SIGNAL(currentItemChanged(QListWidgetItem*, QListWidgetItem*)), this, SLOT(showTableFields(QListWidgetItem*)));

  tableContent = new QTableWidget(10, 2, this);
  QStringList tableLabels;
  tableLabels.append(tr("        Name        "));
  tableLabels.append(tr("        Type        "));
  tableContent->setHorizontalHeaderLabels(tableLabels);
  tableContent->resizeColumnsToContents();
  //tableContent->hide();

  setDB(id);

  QVBoxLayout *layoutListBox = new QVBoxLayout;
  layoutListBox->addWidget(tablesListBox);

  QVBoxLayout *layoutTableContent = new QVBoxLayout;
  layoutTableContent->addWidget(tableContent);

  QHBoxLayout *layout = new QHBoxLayout;
  layout->addLayout(layoutListBox);
  layout->addLayout(layoutTableContent, 1);
  setLayout(layout);
}

void infoTab::destroy_arrays(){
  if (tbl_n_fields){
    for (int i = 0; i < n_tables; i++){
      for (int ii = 0; ii < tbl_n_fields[i]; ii++){
	delete[] tbl_fields_name[i][ii];
	delete[] tbl_fields_type[i][ii];
      }
      delete[] tbl_fields_name[i];
      delete[] tbl_fields_type[i];
    }
    delete[] tbl_fields_name;
    delete[] tbl_fields_type;
    delete[] tbl_n_fields;
    tbl_n_fields = 0;
  }
}

void infoTab::setDB(const int id)
{
  db_id = id;
  destroy_arrays();

  if (tablesListBox){
    sprintf(the_query, "SELECT c.relname as \"Name\", CASE c.relkind WHEN 'r' THEN 'table' WHEN 'v' THEN 'view' WHEN 'i' THEN 'index' WHEN 'S' THEN 'sequence' WHEN 's' THEN 'special' END as \"Type\", r.rolname as \"Owner\" \
                	FROM pg_catalog.pg_class c JOIN pg_catalog.pg_roles r ON r.oid = c.relowner  LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace WHERE c.relkind IN ('r','v','S','') \
	                AND n.nspname NOT IN ('pg_catalog', 'pg_toast') AND pg_catalog.pg_table_is_visible(c.oid) ORDER BY 1,2");
    pg_dbi::get()->query(db_id, the_query);
    pg_dbi::get()->query_dimensions(lines, columns);
    n_tables = lines;

    tablesListBox->clear();
    QStringList tables_list;   //tables_list.clear();
    for (int i = 0; i < n_tables; ++i)
      tables_list.append(tr("%1").arg(pg_dbi::get()->query_value(i,0)));
    tablesListBox->insertItems(0, tables_list);
    for (int i = 0; i < n_tables; ++i) tablesListBox->item(i)->setData(12, i);

    //for (int i = 0; i < n_tables; ++i)
    //  printf("table %s data %d\n", tablesListBox->item(i)->text().toLocal8Bit().constData(), tablesListBox->item(i)->data(12).toInt());

    tbl_n_fields    = new int[n_tables];
    tbl_fields_name = new char**[n_tables];
    tbl_fields_type = new char**[n_tables];

    QString db_tbl_index;
    for (int i = 0; i < n_tables; ++i){
      sprintf(the_query, "SELECT c.oid, n.nspname, c.relname FROM pg_catalog.pg_class c \
                    	  LEFT JOIN pg_catalog.pg_namespace n ON n.oid = c.relnamespace \
			  WHERE c.relname ~ '^(%s)$' AND pg_catalog.pg_table_is_visible(c.oid) ORDER BY 2, 3", tables_list.at(i).toLocal8Bit().constData());
      pg_dbi::get()->query(db_id, the_query);
      db_tbl_index = pg_dbi::get()->query_value(0,0);

      sprintf(the_query, "SELECT a.attname, pg_catalog.format_type(a.atttypid, a.atttypmod), (SELECT substring(pg_catalog.pg_get_expr(d.adbin, d.adrelid) for 128)  FROM pg_catalog.pg_attrdef d WHERE d.adrelid = a.attrelid AND d.adnum = a.attnum AND a.atthasdef), a.attnotnull, a.attnum  \
                 	  FROM pg_catalog.pg_attribute a WHERE a.attrelid = '%d' AND a.attnum > 0 AND NOT a.attisdropped ORDER BY a.attnum;", db_tbl_index.toInt());
      pg_dbi::get()->query(db_id, the_query);
      pg_dbi::get()->query_dimensions(lines, columns);
      tbl_n_fields[i] = lines;

      tbl_fields_name[i] = new char*[tbl_n_fields[i]];
      tbl_fields_type[i] = new char*[tbl_n_fields[i]];
      for (int ii = 0; ii < tbl_n_fields[i]; ++ii){
	tbl_fields_name[i][ii] = new char[256];
	tbl_fields_type[i][ii] = new char[256];
	sprintf(tbl_fields_name[i][ii], pg_dbi::get()->query_value(ii,0));
	sprintf(tbl_fields_type[i][ii], pg_dbi::get()->query_value(ii,1));
	//printf("table %s fname %s\n", tables_list.at(i).toLocal8Bit().constData(), tbl_fields_name[i][ii]);
	//printf("table %s ftype %s\n", tables_list.at(i).toLocal8Bit().constData(), tbl_fields_type[i][ii]);
      }
    }
  }
}

void infoTab::showTableFields(QListWidgetItem *item)
{
  int tbl_id = item->data(12).toInt();
  //printf("database id %d, table id %d\n", db_id, tbl_id);

  if (prnt){
    mainWindow::get()->dispatchMsg(str->sprintf(" >>> TABLE CONTENT: %s\n", tablesListBox->item(tbl_id)->text().toLocal8Bit().constData()), 0);
    mainWindow::get()->dispatchMsg(str->sprintf("                FIELD                         TYPE\n"), 0);
  }
  if (tbl_id >= 0){
    tableContent->clearContents();
    tableContent->setRowCount(tbl_n_fields[tbl_id]);

    for (int i = 0; i < tbl_n_fields[tbl_id]; ++i){
      QTableWidgetItem *tableItemsName = new QTableWidgetItem(tbl_fields_name[tbl_id][i]);
      tableItemsName->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      tableItemsName->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      tableContent->setItem(i, 0, tableItemsName);

      QTableWidgetItem *tableItemsType = new QTableWidgetItem(tbl_fields_type[tbl_id][i]);
      tableItemsType->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
      tableItemsType->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable);
      tableContent->setItem(i, 1, tableItemsType);

      if (prnt)
	mainWindow::get()->dispatchMsg(str->sprintf(" %20s %28s\n", tbl_fields_name[tbl_id][i], tbl_fields_type[tbl_id][i]), 0);
    }
  }

  tableContent->resizeColumnsToContents();
  //tableContent->adjustSize();
  tableContent->update();

  fflush(stdout);
}

queryTab::queryTab(const int id, QWidget *parent)
  : QWidget(parent)
{
  str = new QString;
  if (id > 0) ;
}

