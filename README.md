The package contains functions which can provide sharable output for the models fitted using the gamlss R package. The objects prepared by the functions do not contain any sensitive information protected by the law. (e.g. EGPR).

The list of main functions:

gamlssReport::extracts the paramaters of the fitted GAMLSS object needed to make predictions (uses: extract_terms)

print.gamlssReport::print function for objects created by gamlssReport 

predict.gamlssReport::user friendly wrapper for make_prediction (uses: make_prediction)

centile.gamlssReport:: calculates a centile for y given xs (uses: predict.gamlssReport)

score.gamlssReport:: calculates a score for centile given xs (uses: predict.gamlssReport)

plot.gamlssReport::plots centiles; can also return an object that can be used to make some nicer plots, e.g. using ggplot (uses: score.gamlssReport)



The list of aux functions (not intended to be used by the user):

make_prediction::using the object generated by extract_params, it makes predictions (for the parameters), based on the new data (uses: make_spline)

extract_terms::extract splinevar and fixef formula from the fitted gamlss object

make_spline::makes B-spline basis using splineDesing from splines package. Note: gives the appropriate res only if gamlss was fitted using the default settings.

make_inverse::inverse of the link function. Currently implemented links: log, exp (do we really need that?), logit, inverse, identity