#!/bin/sh
cp `locate createEsetList.R | grep library/curated` . 

R --vanilla "--args patientselection2.config.hold  eset.scaled.rda eset.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2binary.config.hold eset.binary.scaled.rda eset.binary.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2validation.config.hold eset.validation.scaled.rda eset.validation.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2debulking.config.hold eset.debulking.scaled.rda eset.debulking.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2debulkingnoscale.config.hold eset.debulking.rda eset.debulking.log"  < createEsetList.R
R --vanilla "--args patientselection2allos.config.hold eset.allos.scaled.rda eset.allos.scaled.log"  < createEsetList.R
R --vanilla "--args patientselection2inclgrade2.config.hold eset.inclgrade2.scaled.rda eset.inclgrade2.scaled.log"  < createEsetList.R
