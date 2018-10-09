cp photos_gen_analyzer.cc photos_gen_analyzer_muon.cc
find . -name photos_gen_analyzer_muon.cc | xargs perl -pi -e s/photos_gen_analyzer/photos_gen_analyzer_muon/g
find . -name photos_gen_analyzer_muon.cc | xargs perl -pi -e s/"int leppid=11"/"int leppid=13"/g

#int leppid=13
cp photos_gen_analyzer.cc photos_gen_analyzer_electron.cc
find . -name photos_gen_analyzer_electron.cc | xargs perl -pi -e s/photos_gen_analyzer/photos_gen_analyzer_electron/g
find . -name photos_gen_analyzer_electron.cc | xargs perl -pi -e s/"int leppid=13"/"int leppid=11"/g