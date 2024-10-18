#!/usr/bin/env bash

#### ------------------ ####
#### 1. BASE SIMULATION ####
#### ------------------ ####
# 1A. simulation 
Rscript scripts/simulation_model.R \
	--out base_simulation

# 1B. fit base model to base simulation
Rscript scripts/run_stan.R \
	--dat simulations/base_simulation.tsv \
	--stan stan/base_model.stan \
	--seqDesignMatrix "scaled_log10_vl_obs"


#### ------------------- ####
#### 2. FULL  SIMULATION ####
#### ------------------- ####
# 2A. simulate data with FPR and FNR
Rscript scripts/simulation_model.R \
	--prob_MI_fpr 0.01 \
	--prob_MI_fnr 0.3 \
	--out full_simulation

# 2B. fit base model to full simulation
Rscript scripts/run_stan.R \
	--dat simulations/full_simulation.tsv \
	--stan stan/base_model.stan \
	--seqDesignMatrix "scaled_log10_vl_obs"

# 2C. fit full model to full simulation
Rscript scripts/run_stan.R \
	--dat simulations/full_simulation.tsv \
	--stan stan/full_model.stan \
	--seqDesignMatrix "scaled_log10_vl_obs"


#### ------------------------------- ####
#### 3. DELTA SENSITIVITY SIMULATION ####
#### ------------------------------- ####
# simulate across a range of delta bar values (mean prevalence of multiple infection)
for delta in 0 0.05 0.10 0.20; do
	echo $delta
	mkdir -p simulations/sensitivity 
	# 3A. simulate
	sensitivity_simulation=$(Rscript scripts/simulation_model.R \
		--prob_MI_fpr 0.01 \
		--prob_MI_fnr 0.3 \
		--prob_MI_baseline ${delta} \
		--out sensitivity/full_delta${delta})
	# 3B. fit full model to sensitivity simulation
	Rscript scripts/run_stan.R \
		--dat $sensitivity_simulation \
		--stan stan/full_model.stan \
		--seqDesignMatrix "scaled_log10_vl_obs"
done


#### --------------------------------- ####
#### 4. EPSILON SENSITIVITY SIMULATION ####
#### --------------------------------- ####
# simulate across a range of epsilon (false positive rate) values
for fpr in 0 0.005 0.01 0.05; do
	# 4A. simulate
	sensitivity_simulation=$(Rscript scripts/simulation_model.R \
		--prob_MI_fpr $fpr \
		--prob_MI_fnr 0.3 \
		--prob_MI_baseline 0.05 \
		--out sensitivity/full_fpr${fpr})
	# 4B. fit full model to sensitivity simulation
	Rscript scripts/run_stan.R \
		--dat $sensitivity_simulation \
		--stan stan/full_model.stan \
		--seqDesignMatrix "scaled_log10_vl_obs"
done


#### -------------------------------- ####
#### 5. LAMBDA SENSITIVITY SIMULATION ####
#### -------------------------------- ####
# simulate across a range of lambda (false negative rate) values
for fnr in 0.1 0.2 0.3 0.4; do
	# 5A. simulate
	sensitivity_simulation=$(Rscript scripts/simulation_model.R \
		--prob_MI_fpr 0.01 \
		--prob_MI_fnr $fnr \
		--prob_MI_baseline 0.05 \
		--out sensitivity/full_fnr${fnr})
	# 5B. fit full model to sensitivity simulation
	Rscript scripts/run_stan.R \
		--dat $sensitivity_simulation \
		--stan stan/full_model.stan \
		--seqDesignMatrix "scaled_log10_vl_obs"
done


#### ----------------------- ####
#### 6. EXTENDED  SIMULATION ####
#### ----------------------- ####
# simulate a risk factor associated with multiple infection
# 6A. simulate
Rscript scripts/simulation_model.R \
	--prob_MI_fpr 0.01 \
	--prob_MI_fnr 0.3 \
	--null_risk_factors 4 \
	--prob_risk_factor 0.50 \
	--prob_MI_risk_factor 0.10 \
	--out extended_simulation

# 6B. fit extended model to extended simulation with shrinkage priors
Rscript scripts/run_stan.R \
	--dat simulations/extended_simulation.tsv \
	--stan stan/extended_model_hsp.stan \
	--miDesignMatrix \
		risk_factor_MI_1 risk_factor_MI_2 \
		risk_factor_MI_3 risk_factor_MI_4 \
		risk_factor_MI_5 \
	--seqDesignMatrix "scaled_log10_vl_obs" \
	--populationRiskRatios "risk_factor_MI_1=1:risk_factor_MI_1=0"


