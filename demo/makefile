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
PROVENANCE_CMD = $(RSCRIPT_BIN) --vanilla -e 'cat("R.version.string=", R.version.string, "\n", sep=""); blas <- sessionInfo()$$BLAS; if (is.null(blas) || length(blas) == 0L || is.na(blas[1L]) || !nzchar(blas[1L])) blas <- "unknown"; cat("BLAS=", as.character(blas[1L]), "\n", sep="")'

all: run

run: run-$(MODE)

run-serial:
	@for d in $(DEMOS); do \
		echo "[serial] $$d"; \
		out_file=$${d}_serial.Rout; \
		tmp_file=$${out_file}.tmp; \
		if R CMD BATCH --no-save ../$${d}_serial.R "$$tmp_file"; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			: ; \
		else \
			rm -f "$$tmp_file"; \
			echo "[serial] $$d FAILED"; \
			exit 1; \
		fi; \
	done

run-attach:
	@for d in $(DEMOS); do \
		echo "[attach][np=$(NP)] $$d"; \
		out_file=$${d}_npRmpi_attach.Rout; \
		tmp_file=$${out_file}.tmp; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER '' -env R_PROFILE '' $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_N) -n $(NP) $(RSCRIPT_BIN) --no-save ../$${d}_npRmpi_attach.R > "$$tmp_file" 2>&1; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			: ; \
		else \
			rm -f "$$tmp_file"; \
			echo "[attach][np=$(NP)] $$d FAILED"; \
			exit 1; \
		fi; \
	done

run-profile:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(NP)] $$d"; \
		out_file=$${d}_npRmpi_profile.Rout; \
		tmp_file=$${out_file}.tmp; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER $(RPROFILE) -env R_PROFILE '' -env NP_RMPI_PROFILE_RECV_TIMEOUT_SEC $(TIMEOUT_SEC) $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_N) -n $(NP) R CMD BATCH --no-save ../$${d}_npRmpi_profile.R > "$$tmp_file" 2>&1; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			: ; \
		else \
			rm -f "$$tmp_file"; \
			echo "[profile][np=$(NP)] $$d FAILED"; \
			exit 1; \
		fi; \
	done
