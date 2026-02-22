DEMOS = npcdensls npcdensml npcdistls npcmstest npconmode npcopula npdeneqtest npdeptest npindexich npindexks npplreg npqreg npregiv npreglcaic npreglcls npregllaic npregllls npscoef npsdeptest npsigtest npsymtest npudensls npudensml npudistcdf npunitest
RPROFILE ?= $(shell Rscript -e 'cat(system.file(".Rprofile", package="npRmpi"))')

all:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(np)] $$d"; \
		R_PROFILE_USER=$(RPROFILE) mpiexec -n $(np) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R; \
	done
