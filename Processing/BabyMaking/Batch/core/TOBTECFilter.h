
int TOBTEC_ok() {

    int decision = 1;
    for (unsigned int it = 0; it<jets_AK5PF_pt->size(); it++) {
        if (fabs(jets_AK5PF_eta->at(it)) > 0.9 && fabs(jets_AK5PF_eta->at(it)) < 1.9 && jets_AK5PF_chg_Mult->at(it) - jets_AK5PF_neutral_Mult->at(it) > 40) decision = 0;
    }

    return decision;
}

