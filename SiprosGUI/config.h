#ifndef CONFIG_H
#define CONFIG_H
#include <QString>
#include <QList>
#include "mainwindow.h"

class config
{
public:
    config();
    QString searchName;
    QString searchType;
    QStringList fastaDB;
    QStringList rawFiles;
    QString fragmentationMethod;
    QList<int> parentMassWindows;
    int Minimum_Peptide_Length;
    int Maximum_Peptide_Length;
    float Search_Mass_Tolerance_Parent_Ion;
    float Mass_Tolerance_Fragment_Ions;
    QString Cleave_After_Residues;
    QString Cleave_Before_Residues;
    int Maximum_Missed_Cleavages;
    bool Remove_First_Methionine;
    int Max_PTM_Count;
    QString PTM_table;
    QString Isotopic_table;
    QString Training_Decoy_Prefix;
    QString Testing_Decoy_Prefix;
    QString FDR_Filtering_level;
    float FDR_Threshold;
    int Min_Peptide_Per_Protein;
    int Min_Unique_Peptide_Per_Protein;
    float Filter_Mass_Tolerance_Parent_Ion;
    QString Filter_Mass_Tolerance_Parent_Ion_Unit;

    void readJson();
    QString writeJson();
    void readMainWindow(Ui::MainWindow *ui);
    static config& getConfig();
    config(const config& other) = delete;
    config& operator=(const config& other) = delete;
};

#endif // CONFIG_H
