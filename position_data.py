import csv
import re
import sys
import os
import pandas as pd


lanl_features_file = "RefData/LANL_EnvFeatures_Dec7.xlsx"
lanl_features_notes = "RefData/LANL_Features_Dec7_RS_Notes.xlsx"
lanl_features_df = pd.read_excel(lanl_features_file, header=[0])
lanl_features_df = lanl_features_df.rename(columns={"Env feature(s)": "Features"})

lanl_features_notes_df = pd.read_excel(lanl_features_notes, header=[0])
lanl_features_notes_df = lanl_features_notes_df[lanl_features_notes_df["Notes"] == "Y"]
#print(lanl_features_notes_df)





def get_env_features(lanl_features_df, print_dicts = False):
    position_features_df = lanl_features_df.groupby(["Position", "Features"]).size().reset_index(name="Count")
    position_features_df["Features"] = position_features_df["Features"].str.replace(r"gp120, |gp41, ", "", regex=True).str.strip()

    position_to_feature = dict(zip(position_features_df["Position"], position_features_df["Features"]))
    for key, value in list(position_to_feature.items()):
        if value == "gp120" or value == "gp41":
            del position_to_feature[key]      

    feature_to_positions = position_features_df.groupby("Features")["Position"].apply(list).to_dict()
    del feature_to_positions["gp41"]
    del feature_to_positions["gp120"]

    if print_dicts:
        for key, value in position_to_feature.items():
            print(f"{key}: {value}") 
        for key, value in feature_to_positions.items():
            print(f"{key}: {value}")

    return position_to_feature, feature_to_positions

def get_cd4_contacts(lanl_features_df, print_list = False):
    cd4_contacts_df = lanl_features_df[lanl_features_df["Title"] == "CD4 contacts"]
    cd4_contacts = sorted(set(cd4_contacts_df["Position"].tolist()))
    cd4_contacts_refs = sorted(set(cd4_contacts_df["Reference"].tolist()))
    cd4_contacts_methods = sorted(set(cd4_contacts_df["Experimental method(s)"].tolist()))

    if print_list:
        print(cd4_contacts, "\n", cd4_contacts_refs, "\n", cd4_contacts_methods, "\n")

    return cd4_contacts, cd4_contacts_refs, cd4_contacts_methods


def get_bnab_contacts(df, bnab, print_list = False):
    bnab_contacts_df = df[df["Title"].str.contains(bnab) & df["Title"].str.contains("contacts")]
    bnab_contacts = sorted(set(bnab_contacts_df["Position"].tolist()))
    references = sorted(set(bnab_contacts_df["Reference"].tolist()))
    return bnab_contacts, references

def get_escape_positions(df, bnab):
    search_strings = ["immunotherapy", "escape"]
    bnab_escape_df = df[df["Title"].str.contains(bnab) & 
                        df["Title"].str.contains('|'.join(search_strings))]  
    escape_positions = sorted(set(bnab_escape_df["Position"].tolist()))
    references = sorted(set(bnab_escape_df["Reference"].tolist()))
    return escape_positions, references

position_to_feature, feature_to_position = get_env_features(lanl_features_df)
cd4_contacts, cd4_contacts_refs, cd4_contacts_methods = get_cd4_contacts(lanl_features_df)
#print(cd4_contacts, "\n", cd4_contacts_refs, "\n", cd4_contacts_methods, "\n")

bnabs = {
    "cd4_bs": ["VRC01", "3BNC117", "N6", "VRC07-523"],
    "v2": ["PGDM1400", "CAP256"],
    "v3": ["10-1074", "PGT121"],
    "mper": ["10E8"],
    "gp41_interface": ["8ANC195", "PGT151"]
}

for bnab_type, bnab_list in bnabs.items():
    print(f"Processing items in {bnab_type}:")
    for bnab in bnab_list:
        bnab_contacts, contact_refs = get_bnab_contacts(lanl_features_df, bnab)
        escape_positions, escape_refs = get_escape_positions(lanl_features_df, bnab)
        #print(bnab, " contact: ", bnab_contacts, " (", contact_refs, ")")
        print(bnab, " escape: ", escape_positions, " (", escape_refs, ")")

