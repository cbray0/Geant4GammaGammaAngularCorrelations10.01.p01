if [[ -d "$(processes/hadronic/models/)" && -d "$(processes/hadronic/models/de_excitation/)" && -d "$(particles/management)" ]] ; then
    rsync -avh radioactive_decay processes/hadronic/models/
    rsync -avh photon_evaporation processes/hadronic/models/de_excitation/
    rsync -avh management particles/management
else
    echo "Please run from properly configured \"source\" directory. (The folders processes/hadronic/models/, processes/hadronic/models/de_excitation/, and particles/management should exist inside the \"source\" directory.)"
fi
