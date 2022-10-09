The codes is for paper of "STAR-IOS Aided NOMA Networks: Channel Model Approximation and Performance Analysis"
 in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 6861-6876, Sept. 2022, doi: 10.1109/TWC.2022.3152703.

IEEE explore Link: https://ieeexplore.ieee.org/document/9722712

The file named as "Simulation.m" is the outage probability of the simulation results for the three channel models from Lemma 1 to Lemma 3.

The file named as "Analysis_center_limit_model.m" is the analytical outage pribability under the center limit model (Theorems 1-2)

The file named as "Analysis_curve_fitting.m"  is the analytical outage pribability under the curve fitting model (Theorems 3-4)

The file named as "Analysis_N_fold_as_diversity_orders.m" is the analytical outage pribability under the M-fold convolution model (Theorems 5-6, and Corollaries 2-7).

The file named as "Curve_fitting_coefficients.m" is to present how we obtain the coefficients in the curve fitting model (Lemma 3)

---------------------------------------------------------------------------
Please first run the simulation results for three models:
./Simulation.m

For validation for the the center limit model, please use:
./Analysis_center_limit_model.m

For validation for the curve fitting model, please use:
1) ./Curve_fitting_coefficients.m  --> this is for the coefficients of the  
2) ./Analysis_curve_fitting.m

For validation for the M-fold convolution model and to investigate the diversity orders, please use:
./Analysis_N_fold_as_diversity_orders.m
