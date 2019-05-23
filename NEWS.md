# intrinsicDimension 1.2.0
## New features
- The reference values for the ESS dimension estimation can now be accessed by the function `essReference`.
## Minor improvements and fixes
- The reference values for ESS dimension estimation is computed with an algorithm that avoids gamma function evaluations and is thus numerically stable for much higher dimensions.