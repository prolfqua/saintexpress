# Changelog

## saintexpress 0.99.0

- Licence changed from MIT to GPL (\>= 3). The package is a derivative
  work of the GPL-3-or-later SAINTexpress C++ sources, so a permissive
  licence was not applicable. Hyungwon Choi and Damian Fermin are now
  credited as authors and copyright holders.

- Attribution added: the README now documents that this package is an
  AI-assisted port of the GPL-3 licensed SAINTexpress C++ sources by
  Hyungwon Choi and Damian Fermin.

- BiocStyle moved from Imports to Suggests; it is only used to build the
  vignette.

- Author credits updated: Alexey Nesvizhskii is credited as an author,
  Witold Wolski is now the maintainer.

- Intensity scoring now handles constant complete control profiles with
  the median-control-SD fallback used by native SAINTexpress, while
  reporting the affected prey identifiers.

- Initial pure-R implementation of SAINTexpress spectral-count and
  intensity scoring with SAINT input validation.
