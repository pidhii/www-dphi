.PHONY: MC data MC_incl

LISTS_MC := \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Ariadne_Low_Q2_NC_DIS_040506e.list \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Ariadne_Low_Q2_NC_DIS_0607p.list

LISTS_MC_INCL := \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_ari_incl_nc_DIS_lowQ2_040506e.list \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_ari_incl_nc_DIS_lowQ2_0607p.list

LISTS_LEPTO := \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Lepto_low_Q2_NC_DIS_040506e.list \
		/afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_Lepto_low_Q2_NC_DIS_0607p.list

LISTS_DATA := \
		$(foreach list,05e 06e 0607p, /afs/desy.de/user/q/quintera/public/Jets/Lists/Sample_$(list).list)

ifeq ($(out),)
$(error Unspecified output-directory. (can set with `out=...`))
endif

ifdef jet_minEt
export CUT_JET_MINET=$(jet_minEt)
$(warning --- Jet min Et = $(CUT_JET_MINET))
endif

ifdef lepton_minE
export CUT_LEPTON_MINE=$(lepton_minE)
$(warning --- Lepton min E = $(CUT_LEPTON_MINE))
endif

ifdef lepton_miny
export CUT_LEPTON_MINY=$(lepton_miny)
$(warning --- Lepton min y = $(CUT_LEPTON_MINY))
endif

ifdef lepton_maxy
export CUT_LEPTON_MAXY=$(lepton_maxy)
$(warning --- Lepton max y = $(CUT_LEPTON_MAXY))
endif

ifdef reweight
#
# Check if sources are set.
#
ifndef data
$(error Data source for reweighting unsepcified)
endif
ifndef reco
$(error Reco source for reweighting unsepcified)
endif

export REWEIGHT=yes
export REWEIGHT_DATA=$(data)
export REWEIGHT_RECO=$(reco)

$(warning --- Reweighting enabled.)
endif

ifdef macro
MACRO := $(macro)
else
MACRO := MakeHist.C
$(warning Using default ROOT macro - "$(MACRO)")
endif

MC:
	for list in $(LISTS_MC); do ./bin/run.pl --macro $(MACRO) -o $(out) --list $$list; done

MC_incl:
	for list in $(LISTS_MC_INCL); do ./bin/run.pl --macro $(MACRO) -o $(out) --list $$list; done

lepto:
	for list in $(LISTS_LEPTO); do ./bin/run.pl --macro $(MACRO) -o $(out) --list $$list; done

data:
	for list in $(LISTS_DATA); do ./bin/run.pl --macro $(MACRO) -o $(out) --list $$list; done
