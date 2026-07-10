DEMOS ?= npcdens npcdist npcmstest npconmode npcopula npdeneqtest npdeptest npindex npplreg npqreg npreg npregiv npscoef npsdeptest npsigtest npsymtest npudens npudist npunitest

MODE ?= attach
NP ?= 2
NSLAVES ?= 1
TIMEOUT_SEC ?= 600
CONTINUE_ON_ERROR ?= false
DEMO_SRC ?= ..
NP_DEMO_TIER ?= smoke
NP_DEMO_CASES ?=
NP_DEMO_MATRIX ?=
RPROFILE ?= $(if $(wildcard $(DEMO_SRC)/../inst/Rprofile),$(DEMO_SRC)/../inst/Rprofile,$(if $(wildcard $(DEMO_SRC)/../.Rprofile),$(DEMO_SRC)/../.Rprofile,$(shell Rscript -e 'cat(system.file("Rprofile", package="npRmpi"))')))
RSCRIPT_BIN ?= $(shell R RHOME)/bin/Rscript
SETSID_BIN ?= $(shell command -v setsid 2>/dev/null)
SESSION_LAUNCH = $(if $(SETSID_BIN),$(SETSID_BIN),)
MPI_ENV_NP_DEMO_N = $(if $(NP_DEMO_N),-env NP_DEMO_N $(NP_DEMO_N),)
MPI_ENV_NP_DEMO_N_FRAC = $(if $(NP_DEMO_N_FRAC),-env NP_DEMO_N_FRAC $(NP_DEMO_N_FRAC),)
MPI_ENV_NP_DEMO_SRC = -env NP_DEMO_SRC $(DEMO_SRC)
MPI_ENV_NP_DEMO_TIER = -env NP_DEMO_TIER $(NP_DEMO_TIER)
MPI_ENV_NP_DEMO_CASES = $(if $(NP_DEMO_CASES),-env NP_DEMO_CASES "$(NP_DEMO_CASES)",)
MPI_ENV_NP_DEMO_MATRIX = $(if $(NP_DEMO_MATRIX),-env NP_DEMO_MATRIX $(NP_DEMO_MATRIX),)
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
		if NP_DEMO_SRC=$(DEMO_SRC) NP_DEMO_TIER=$(NP_DEMO_TIER) NP_DEMO_CASES="$(NP_DEMO_CASES)" NP_DEMO_MATRIX="$(NP_DEMO_MATRIX)" R CMD BATCH --no-save $(DEMO_SRC)/$${d}_serial.R "$$tmp_file"; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			: ; \
		else \
			status=$$?; \
			{ $(PROVENANCE_CMD); test -s "$$tmp_file" && cat "$$tmp_file"; echo "DEMO_STATUS demo=$$d mode=serial status=failed exit_code=$$status timeout_sec=NA"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			echo "[serial] $$d FAILED"; \
			if test "$(CONTINUE_ON_ERROR)" = "true"; then continue; else exit 1; fi; \
		fi; \
	done

run-attach:
	@for d in $(DEMOS); do \
		echo "[attach][np=$(NP)] $$d"; \
		out_file=$${d}_npRmpi_attach.Rout; \
		tmp_file=$${out_file}.tmp; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER '' -env R_PROFILE '' $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_SRC) $(MPI_ENV_NP_DEMO_TIER) $(MPI_ENV_NP_DEMO_CASES) $(MPI_ENV_NP_DEMO_MATRIX) $(MPI_ENV_NP_DEMO_N) $(MPI_ENV_NP_DEMO_N_FRAC) -n $(NP) $(RSCRIPT_BIN) --no-save $(DEMO_SRC)/$${d}_npRmpi_attach.R > "$$tmp_file" 2>&1; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			: ; \
		else \
			status=$$?; \
			{ $(PROVENANCE_CMD); test -s "$$tmp_file" && cat "$$tmp_file"; echo "DEMO_STATUS demo=$$d mode=attach np=$(NP) status=failed exit_code=$$status timeout_sec=$(TIMEOUT_SEC)"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			echo "[attach][np=$(NP)] $$d FAILED"; \
			if test "$(CONTINUE_ON_ERROR)" = "true"; then continue; else exit 1; fi; \
		fi; \
	done

run-profile:
	@for d in $(DEMOS); do \
		echo "[profile][np=$(NP)] $$d"; \
		out_file=$${d}_npRmpi_profile.Rout; \
		batch_file=$${out_file}.batch; \
		tmp_file=$${out_file}.tmp; \
		if /opt/local/bin/timeout $(TIMEOUT_SEC) mpiexec -env R_PROFILE_USER $(RPROFILE) -env R_PROFILE '' -env NP_RMPI_PROFILE_RECV_TIMEOUT_SEC $(TIMEOUT_SEC) $(MPI_ENV_FI) $(MPI_ENV_NP_DEMO_SRC) $(MPI_ENV_NP_DEMO_TIER) $(MPI_ENV_NP_DEMO_CASES) $(MPI_ENV_NP_DEMO_MATRIX) $(MPI_ENV_NP_DEMO_N) $(MPI_ENV_NP_DEMO_N_FRAC) -n $(NP) R CMD BATCH --no-save $(DEMO_SRC)/$${d}_npRmpi_profile.R "$$batch_file" > "$$tmp_file" 2>&1; then \
			if test ! -s "$$batch_file"; then \
				{ $(PROVENANCE_CMD); test -s "$$tmp_file" && cat "$$tmp_file"; echo "DEMO_STATUS demo=$$d mode=profile np=$(NP) status=failed exit_code=empty_batch timeout_sec=$(TIMEOUT_SEC)"; } > "$$out_file"; \
				rm -f "$$tmp_file" "$$batch_file"; \
				echo "[profile][np=$(NP)] $$d FAILED (empty batch transcript)"; \
				if test "$(CONTINUE_ON_ERROR)" = "true"; then continue; else exit 1; fi; \
			fi; \
			{ $(PROVENANCE_CMD); cat "$$batch_file"; if test -s "$$tmp_file"; then cat "$$tmp_file"; fi; } > "$$out_file"; \
			rm -f "$$tmp_file" "$$batch_file"; \
			: ; \
		else \
			status=$$?; \
			{ $(PROVENANCE_CMD); test -s "$$batch_file" && cat "$$batch_file"; test -s "$$tmp_file" && cat "$$tmp_file"; echo "DEMO_STATUS demo=$$d mode=profile np=$(NP) status=failed exit_code=$$status timeout_sec=$(TIMEOUT_SEC)"; } > "$$out_file"; \
			rm -f "$$tmp_file" "$$batch_file"; \
			echo "[profile][np=$(NP)] $$d FAILED"; \
			if test "$(CONTINUE_ON_ERROR)" = "true"; then continue; else exit 1; fi; \
		fi; \
	done

run-session:
	@for d in $(DEMOS); do \
		echo "[session][nslaves=$(NSLAVES)] $$d"; \
		out_file=$${d}_npRmpi_session.Rout; \
		tmp_file=$${out_file}.tmp; \
		if ! test -f $(DEMO_SRC)/$${d}_npRmpi_session.R; then \
			echo "[session][nslaves=$(NSLAVES)] $$d SKIPPED (no session demo script)"; \
			continue; \
		fi; \
		if $(SESSION_LAUNCH) env NP_DEMO_SRC=$(DEMO_SRC) NP_DEMO_TIER=$(NP_DEMO_TIER) NP_DEMO_CASES="$(NP_DEMO_CASES)" NP_DEMO_MATRIX="$(NP_DEMO_MATRIX)" NP_DEMO_N=$(NP_DEMO_N) NP_DEMO_N_FRAC=$(NP_DEMO_N_FRAC) NP_DEMO_NSLAVES=$(NSLAVES) $(RSCRIPT_BIN) --vanilla $(DEMO_SRC)/$${d}_npRmpi_session.R > "$$tmp_file" 2>&1; then \
			{ $(PROVENANCE_CMD); cat "$$tmp_file"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			pgrep -f 'slavedaemon\.R|Rslaves\.sh' | xargs kill >/dev/null 2>&1 || true; \
			: ; \
		else \
			status=$$?; \
			{ $(PROVENANCE_CMD); test -s "$$tmp_file" && cat "$$tmp_file"; echo "DEMO_STATUS demo=$$d mode=session nslaves=$(NSLAVES) status=failed exit_code=$$status timeout_sec=NA"; } > "$$out_file"; \
			rm -f "$$tmp_file"; \
			pgrep -f 'slavedaemon\.R|Rslaves\.sh' | xargs kill >/dev/null 2>&1 || true; \
			echo "[session][nslaves=$(NSLAVES)] $$d FAILED"; \
			if test "$(CONTINUE_ON_ERROR)" = "true"; then continue; else exit 1; fi; \
		fi; \
	done
