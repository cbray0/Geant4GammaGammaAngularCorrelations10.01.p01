#!/bin/bash
if [[ -d "processes/hadronic/models/" && -d "processes/hadronic/models/de_excitation/" && -d "particles/management" ]] ; then
    rsync -avh ~/Geant4GammaGammaAngularCorrelations10.02.p03/radioactive_decay processes/hadronic/models/
    rsync -avh ~/Geant4GammaGammaAngularCorrelations10.02.p03/photon_evaporation processes/hadronic/models/de_excitation/
else
    echo "Please run from properly configured \"source\" directory. (The folders processes/hadronic/models/, processes/hadronic/models/de_excitation/, and particles/management should exist inside the \"source\" directory.)"
fi
