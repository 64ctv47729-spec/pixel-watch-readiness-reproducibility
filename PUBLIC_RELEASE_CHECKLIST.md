# Public release checklist

The package is technically prepared as a public-release candidate, but it should not be uploaded until each item below is confirmed.

- [x] The participant information sheet and consent wording permit public sharing of pseudonymized repeated health data.
- [x] The ethics approval and institutional data-governance requirements permit this release.
- [x] A disclosure-risk review has considered the combination of repeated wellness, physiology and device scores.
- [x] All authors have approved the package contents and creator order.
- [x] A code licence and, if data are released, a data licence have been selected with institutional approval.
- [x] The repository name, description and final citation metadata have been agreed.
- [ ] The final manuscript and supporting tables have been compared with `tables/` after all author edits.
- [x] `./validate_package.sh` completes successfully.
- [x] The SHA-256 manifest has been regenerated after the last permitted change.
- [ ] A versioned release has been created and archived in a DOI-minting repository.
- [ ] The manuscript Data Availability Statement cites the final DOI or persistent URL and accurately describes any restrictions.

If public release of the analytic data is not permitted, publish the code, documentation, tables and aggregate outputs without `data/analytic_dataset.csv`; explain the restriction in the Data Availability Statement and use an approved controlled-access process where appropriate.

