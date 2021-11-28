#include "config.h"
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QButtonGroup>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>

config::config()
{

}

config& config::getConfig(){
    static config mConfig;
    return mConfig;
}

void config::readMainWindow(Ui::MainWindow *ui)
{
    searchName=ui->lineEditSearchName->text().trimmed();

    QButtonGroup *buttonGroupSeachType=new QButtonGroup;
    buttonGroupSeachType->addButton(ui->radioButtonDDA,1);
    buttonGroupSeachType->addButton(ui->radioButtonDIA,2);
    buttonGroupSeachType->addButton(ui->radioButtonSIP,3);
    buttonGroupSeachType->addButton(ui->radioButtonMutation,4);
    switch (buttonGroupSeachType->checkedId())
    {
        case 1: searchType="DDA";break;
        case 2: searchType="DIA";break;
        case 3: searchType="SIP";break;
        case 4: searchType="Mutation";break;
    };

    fastaDB.clear();
    QTableWidget *fasta=ui->tableWidgetFasta;
    int nFasta=fasta->rowCount();
    for (int i=0;i<nFasta;i++)
    {
        fastaDB.append(fasta->model()->index(i,0).data().toString());
    }

    rawFiles.clear();
    QTableWidget *raw=ui->tableWidgetRaw;
    int nRaw=raw->rowCount();
    for (int i=0;i<nRaw;i++)
    {
        rawFiles.append(raw->model()->index(i,0).data().toString());
    }

    QButtonGroup *buttonGroupFragmentationMethod=new QButtonGroup;
    buttonGroupFragmentationMethod->addButton(ui->radioButtonCID,1);
    buttonGroupFragmentationMethod->addButton(ui->radioButtonHCD,2);
    buttonGroupFragmentationMethod->addButton(ui->radioButtonETD,3);
    switch (buttonGroupFragmentationMethod->checkedId())
    {
        case 1: fragmentationMethod="CID";break;
        case 2: fragmentationMethod="HCD";break;
        case 3: fragmentationMethod="ETD";break;
    };

    parentMassWindows.clear();
    if(ui->checkBoxPminus1->isChecked()) parentMassWindows.append(-1);
    if(ui->checkBoxP0->isChecked()) parentMassWindows.append(0);
    if(ui->checkBoxP1->isChecked()) parentMassWindows.append(1);
    if(ui->checkBoxP2->isChecked()) parentMassWindows.append(2);
    if(ui->checkBoxP3->isChecked()) parentMassWindows.append(3);

    Minimum_Peptide_Length=ui->spinBoxPepLenMin->value();
    Maximum_Peptide_Length=ui->spinBoxPepLenMax->value();

    Search_Mass_Tolerance_Parent_Ion=ui->doubleSpinBoxPrecursorTolerance->value();
    Mass_Tolerance_Fragment_Ions=ui->doubleSpinBoxFragmentTolerance->value();

    Cleave_After_Residues=ui->lineEditCleaveAfterResidue->text().trimmed();
    Cleave_Before_Residues=ui->lineEditCleaveBeforeResidue->text().trimmed();

    QButtonGroup *buttonGroupRemove_First_Methionine=new QButtonGroup;
    buttonGroupRemove_First_Methionine->addButton(ui->radioButtonRemoveFirstMno,1);
    buttonGroupRemove_First_Methionine->addButton(ui->radioButtonRemoveFirstMyes,2);
    switch (buttonGroupRemove_First_Methionine->checkedId())
    {
        case 1: Remove_First_Methionine=false;break;
        case 2: Remove_First_Methionine=true;break;
    };

    Max_PTM_Count=ui->spinBoxPTMmaxNumber->value();

    Testing_Decoy_Prefix=ui->lineEditTestDecoy->text().trimmed();
    Training_Decoy_Prefix=ui->lineEditTarinDecoy->text().trimmed();

    QButtonGroup *buttonGroupFDR_Filtering_level=new QButtonGroup;
    buttonGroupFDR_Filtering_level->addButton(ui->radioButtonFDRpepetide,1);
    buttonGroupFDR_Filtering_level->addButton(ui->radioButtonFDRPSM,2);
    switch (buttonGroupFDR_Filtering_level->checkedId())
    {
        case 1: FDR_Filtering_level="Peptide";break;
        case 2: FDR_Filtering_level="PSM";break;
    };

    FDR_Threshold=ui->doubleSpinBoxFDRthreshold->value();

    Min_Peptide_Per_Protein=ui->spinBoxMinPeptide->value();

    Min_Unique_Peptide_Per_Protein=ui->spinBoxMinUniquePetide->value();

    Filter_Mass_Tolerance_Parent_Ion=ui->doubleSpinBoxPrecursorMassTolerance->value();

    QButtonGroup *buttonGroupFilter_Parent_Ion_Unit=new QButtonGroup;
    buttonGroupFilter_Parent_Ion_Unit->addButton(ui->radioButtonMassErrorPPM,1);
    buttonGroupFilter_Parent_Ion_Unit->addButton(ui->radioButtonMassErrorDa,2);
    switch (buttonGroupFilter_Parent_Ion_Unit->checkedId())
    {
        case 1: Filter_Mass_Tolerance_Parent_Ion_Unit="PPM";break;
        case 2: Filter_Mass_Tolerance_Parent_Ion_Unit="Da";break;
    };
}

QString config::writeJson(){
    QJsonObject json;
    QJsonArray jsonArray;
    QJsonDocument jsonDoc;
    QString doc;

    json.insert("searchName",searchName);
    json.insert("searchType",searchType);

    for (QString db: fastaDB)
    {
        jsonArray.append(db);
    }
    json.insert("fastaDB",jsonArray);

    jsonArray=QJsonArray();
    for (QString db: rawFiles)
    {
        jsonArray.append(db);
    }
    json.insert("rawFiles",jsonArray);

    json.insert("fragmentationMethod",fragmentationMethod);

    jsonArray=QJsonArray();
    for (int window: parentMassWindows)
    {
        jsonArray.append(window);
    }
    json.insert("parentMassWindows",jsonArray);

    json.insert("Minimum_Peptide_Length",Minimum_Peptide_Length);
    json.insert("Maximum_Peptide_Length",Maximum_Peptide_Length);
    json.insert("Search_Mass_Tolerance_Parent_Ion",Search_Mass_Tolerance_Parent_Ion);
    json.insert("Search_Mass_Tolerance_Parent_Ion",Search_Mass_Tolerance_Parent_Ion);
    json.insert("Mass_Tolerance_Fragment_Ions",Mass_Tolerance_Fragment_Ions);
    json.insert("Cleave_After_Residues",Cleave_After_Residues);
    json.insert("Cleave_Before_Residues",Cleave_Before_Residues);
    json.insert("Maximum_Missed_Cleavages",Maximum_Missed_Cleavages);
    json.insert("Remove_First_Methionine",Remove_First_Methionine);
    json.insert("Max_PTM_Count",Max_PTM_Count);
    json.insert("PTM_table",PTM_table);
    json.insert("Isotopic_table",Isotopic_table);
    json.insert("Training_Decoy_Prefix",Training_Decoy_Prefix);
    json.insert("Testing_Decoy_Prefix",Testing_Decoy_Prefix);
    json.insert("FDR_Filtering_level",FDR_Filtering_level);
    json.insert("FDR_Threshold",FDR_Threshold);
    json.insert("Min_Peptide_Per_Protein",Min_Peptide_Per_Protein);
    json.insert("Min_Unique_Peptide_Per_Protein",Min_Unique_Peptide_Per_Protein);
    json.insert("Filter_Mass_Tolerance_Parent_Ion",Filter_Mass_Tolerance_Parent_Ion);
    json.insert("Filter_Mass_Tolerance_Parent_Ion_Unit",Filter_Mass_Tolerance_Parent_Ion_Unit);

    jsonDoc.setObject(json);
    doc=jsonDoc.toJson();
    return doc;
}





