# Analysis data

`analytic_dataset.csv` contains the minimum variables needed to reproduce the reported models and sensitivities.

Disclosure reduction applied before packaging:

- direct response identifiers and email-derived crosswalks were removed;
- participant identifiers were remapped sequentially from 1 to 38;
- exact questionnaire, Kubios and Fitabase timestamps were removed;
- calendar dates were replaced with synthetic dates that preserve within-participant ordering, protocol day and episode separation;
- free-text and unused questionnaire fields were excluded;
- the duplicate-adjudication file was reduced to the pseudonymous participant-date keys required by the sensitivity analysis.

The synthetic `q_date` values must not be interpreted as actual collection dates. They are retained solely because the prediction workflow groups repeated predictions by participant-day.

`audit_summary.json` contains aggregate counts used when recreating summary tables. It excludes raw-source filenames, email addresses and real date ranges.

`DATA_DICTIONARY.csv` describes every released variable and its analytic role.

