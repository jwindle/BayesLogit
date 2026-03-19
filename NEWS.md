# BayesLogit 2.3 (2026-03-19)

- Updated CITATION to use `bibentry()` replacing deprecated `citEntry()` and `personList()`.
- Added published JASA article to citations.
- Fixed stale and broken URLs (`http://` to `https://`, removed broken thesis repo link).
- Renamed internal functions with dots in names to avoid false S3 method detection.
- Updated package metadata (version, date, URL field).

# BayesLogits 2.2 (2025-06-27)

-Maintenance release: updated CITATION to use `bibentry()`, fixed stale URLs,
-added published JASA article to citations, and updated package metadata.

# BayesLogit 2.1 (2019-09-25)

- BayesLogit has been completely re-organized. Previously, the package was
  derived from a repository used for a doctoral thesis, and included many
  items beyond what was necessary for an R package. This made the
  subsequent package difficult to maintain.
  The new package and repository (<https://github.com/jwindle/BayesLogit>)
  has been slimmed down and now is focused on implementation of PolyaGamma
  sampling methods.
  Thanks go to Jari Oksanen and James Balamuta for their help in getting
  thepackage ready for CRAN. 
