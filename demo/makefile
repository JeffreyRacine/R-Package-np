DEMOS = npcdensls npcdensml npcdistls npcmstest npconmode npcopula npdeneqtest npdeptest npindexich npindexks npplreg npqreg npregiv npreglcaic npreglcls npregllaic npregllls npscoef npsdeptest npsigtest npsymtest npudensls npudensml npudistcdf npunitest

MODE ?= attach
NP ?= 2
RPROFILE ?= $(if $(wildcard ../.Rprofile),../.Rprofile,$(shell Rscript -e 'cat(system.file("Rprofile", package="npRmpi"))'))
MPI_ENV_NP_DEMO_N = $(if $(NP_DEMO_N),-env NP_DEMO_N $(NP_DEMO_N),)

all: run

run: run-$(MODE)

run-serial:
	@for d in $(DEMOS); do \
		echo "[serial] $$d"; \
		R CMD BATCH --vanilla ../$${d}_serial.R; \
	done

run-attach:
	@for d in $(DEMOS); do \
		echo "[attach][np=$(NP)] $$d"; \
		mpiexec $(MPI_ENV_NP_DEMO_N) -n $(NP) Rscript ../$${d}_npRmpi_attach.R > $${d}_npRmpi_attach.Rout; \
	done

run-profile:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(NP)] $$d"; \
		R_PROFILE_USER=$(RPROFILE) mpiexec $(MPI_ENV_NP_DEMO_N) -n $(NP) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R; \
	done
