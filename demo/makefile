DEMOS ?= npcdensls npcdensml npcdistls npcmstest npconmode npcopula npdeneqtest npdeptest npindexich npindexks npplreg npqreg npregiv npreglcaic npreglcls npregllaic npregllls npscoef npsdeptest npsigtest npsymtest npudensls npudensml npudistcdf npunitest

MODE ?= attach
NP ?= 2
TIMEOUT_SEC ?= 180
RPROFILE ?= $(if $(wildcard ../inst/Rprofile),../inst/Rprofile,$(if $(wildcard ../.Rprofile),../.Rprofile,$(shell Rscript -e 'cat(system.file("Rprofile", package="npRmpi"))')))
RSCRIPT_BIN ?= $(shell R RHOME)/bin/Rscript
MPI_ENV_NP_DEMO_N = $(if $(NP_DEMO_N),-env NP_DEMO_N $(NP_DEMO_N),)

all: run

run: run-$(MODE)

run-serial:
	@for d in $(DEMOS); do \
		echo "[serial] $$d"; \
		R CMD BATCH --no-save ../$${d}_serial.R; \
	done

run-attach:
	@for d in $(DEMOS); do \
		echo "[attach][np=$(NP)] $$d"; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER '' -env R_PROFILE '' -env FI_TCP_IFACE en0 -env FI_PROVIDER tcp -env FI_SOCKETS_IFACE en0 $(MPI_ENV_NP_DEMO_N) -n $(NP) $(RSCRIPT_BIN) --no-save ../$${d}_npRmpi_attach.R > $${d}_npRmpi_attach.Rout 2>&1; then \
			: ; \
		else \
			echo "[attach][np=$(NP)] $$d retry on lo0"; \
			/opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER '' -env R_PROFILE '' -env FI_TCP_IFACE lo0 -env FI_PROVIDER tcp -env FI_SOCKETS_IFACE lo0 $(MPI_ENV_NP_DEMO_N) -n $(NP) $(RSCRIPT_BIN) --no-save ../$${d}_npRmpi_attach.R > $${d}_npRmpi_attach.Rout 2>&1; \
		fi; \
	done

run-profile:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(NP)] $$d"; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER $(RPROFILE) -env R_PROFILE '' -env NP_RMPI_PROFILE_RECV_TIMEOUT_SEC $(TIMEOUT_SEC) -env FI_TCP_IFACE en0 -env FI_PROVIDER tcp -env FI_SOCKETS_IFACE en0 $(MPI_ENV_NP_DEMO_N) -n $(NP) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R > $${d}_npRmpi_profile.Rout 2>&1; then \
			: ; \
		else \
			echo "[profile][np=$(NP)] $$d retry on lo0"; \
			/opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER $(RPROFILE) -env R_PROFILE '' -env NP_RMPI_PROFILE_RECV_TIMEOUT_SEC $(TIMEOUT_SEC) -env FI_TCP_IFACE lo0 -env FI_PROVIDER tcp -env FI_SOCKETS_IFACE lo0 $(MPI_ENV_NP_DEMO_N) -n $(NP) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R > $${d}_npRmpi_profile.Rout 2>&1; \
		fi; \
	done
