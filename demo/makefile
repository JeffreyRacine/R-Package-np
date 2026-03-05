DEMOS ?= npcdensls npcdensml npcdistls npcmstest npconmode npcopula npdeneqtest npdeptest npindexich npindexks npplreg npqreg npregiv npreglcaic npreglcls npregllaic npregllls npscoef npsdeptest npsigtest npsymtest npudensls npudensml npudistcdf npunitest

MODE ?= attach
NP ?= 2
TIMEOUT_SEC ?= 180
RPROFILE ?= $(if $(wildcard ../inst/Rprofile),../inst/Rprofile,$(if $(wildcard ../.Rprofile),../.Rprofile,$(shell Rscript -e 'cat(system.file("Rprofile", package="npRmpi"))')))
RSCRIPT_BIN ?= $(shell R RHOME)/bin/Rscript
MPI_ENV_NP_DEMO_N = $(if $(NP_DEMO_N),-env NP_DEMO_N $(NP_DEMO_N),)
FI_PROVIDER ?=
FI_TCP_IFACE ?=
FI_SOCKETS_IFACE ?=
MPI_ENV_FI = $(if $(FI_PROVIDER),-env FI_PROVIDER $(FI_PROVIDER),) $(if $(FI_TCP_IFACE),-env FI_TCP_IFACE $(FI_TCP_IFACE),) $(if $(FI_SOCKETS_IFACE),-env FI_SOCKETS_IFACE $(FI_SOCKETS_IFACE),)

all: run

run: run-$(MODE)

run-serial:
	@for d in $(DEMOS); do \
		echo "[serial] $$d"; \
		if R CMD BATCH --no-save ../$${d}_serial.R; then \
			: ; \
		else \
			echo "[serial] $$d FAILED"; \
			exit 1; \
		fi; \
	done

run-attach:
	@for d in $(DEMOS); do \
		echo "[attach][np=$(NP)] $$d"; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER '' -env R_PROFILE '' $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_N) -n $(NP) $(RSCRIPT_BIN) --no-save ../$${d}_npRmpi_attach.R > $${d}_npRmpi_attach.Rout 2>&1; then \
			: ; \
		else \
			echo "[attach][np=$(NP)] $$d FAILED"; \
			exit 1; \
		fi; \
	done

run-profile:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(NP)] $$d"; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER $(RPROFILE) -env R_PROFILE '' -env NP_RMPI_PROFILE_RECV_TIMEOUT_SEC $(TIMEOUT_SEC) $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_N) -n $(NP) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R > $${d}_npRmpi_profile.Rout 2>&1; then \
			: ; \
		else \
			echo "[profile][np=$(NP)] $$d FAILED"; \
			exit 1; \
		fi; \
	done
