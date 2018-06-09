#ifndef _BX_DB_MON_MAIN_WINDOW_H
#define _BX_DB_MON_MAIN_WINDOW_H

#include "dbMon.h"
#include <QMainWindow>

//QT_BEGIN_NAMESPACE
class QAction;
class QString;
class QLabel;
class QTextBrowser;
class DBDialog;
class BCMDialog;
class LMSDialog;
class MCMDialog;
class MMSDialog;
class DMPDialog;
class DMPMuDialog;
class HVDialog;
// class CHSDialog;
// class CBSDialog;
// class CPSDialog;
// class TMSDialog;
class QSocketNotifier;
//QT_END_NAMESPACE

class mainWindow : public QMainWindow
{
  Q_OBJECT

  public:
    mainWindow();
    ~mainWindow();

    static mainWindow* get() { if (!me) me = new mainWindow; return me; }
    static int current_run() { return cur_run; }

  public slots:
    void dispatchMsg(const QString&, int);

  protected:
    void closeEvent(QCloseEvent *event);

  private slots:
    void openDB();
    void closeDB();
    bool operateDB();

    void bchMon();
    void lmsMon();
    void mchMon();
    void mmsMon();
    void dmpMon();
    void dmmMon();
    void hvMon();
//     void chsMon();
//     void cbsMon();
//     void cpsMon();
//     void tmsMon();

    void help();
    void about();

    void updateStatusBar();
    void monitorModified();

  private:
    pg_dbi *db_mng;
    static mainWindow *me;
    static int cur_run;

    void createActions();
    void createMenus();
    void createContextMenu();
    void createToolBars();
    void createStatusBar();

    QSocketNotifier *stdoutNotif, *stderrNotif;

    QTextBrowser *textBrowser;
    DBDialog  *dbMngDialog;
    BCMDialog *bchMonDialog;
    LMSDialog *lmsMonDialog;
    MCMDialog *mchMonDialog;
    MMSDialog *mmsMonDialog;
    DMPDialog *dmpMonDialog;
    DMPMuDialog *dmpmuMonDialog;
    HVDialog  *hvMonDialog;
//     CHSDialog *chsMonDialog;
//     CBSDialog *cbsMonDialog;
//     CPSDialog *cpsMonDialog;
//     TMSDialog *tmsMonDialog;

    QLabel *statusLabel;
    QLabel *messageLabel;

    QAction *openAllDBAction;
    QAction *openDBActions[NDB];
    QAction *closeDBActions[NDB];
    QAction *operateDBActions[NDB];
    QAction *closeAllDBAction;

    // monitor actions
    QAction *bchMonAction;
    QAction *lmsMonAction;
    QAction *mchMonAction;
    QAction *mmsMonAction;
    QAction *dmpMonAction;
    QAction *dmpmuMonAction;
    QAction *hvMonAction;
//     QAction *chsMonAction;
//     QAction *cbsMonAction;
//     QAction *cpsMonAction;
//     QAction *tmsMonAction;

    QAction *helpAction;
    QAction *aboutAction;
    QAction *aboutQtAction;

    QMenu *dbMenu;
    QMenu *dbButtons[NDB];
    QMenu *monitorMenu;
    QMenu *dumpMenu;
    QMenu *helpMenu;

    //QToolBar *browserToolBar;
};

#endif
