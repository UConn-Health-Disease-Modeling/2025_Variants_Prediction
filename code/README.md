# Variant filtering: 02 to 03

Run the scripts from the project root:

```bash
Rscript code/02_group_variants.R
Rscript code/03_filter_7day.R
```

## Counts after each step

The counts below were calculated from the currently saved 02 outputs and
`variants_who2.rds`.

| Step | Filtering result | Daily rows | Unique Pango variants | Country–variant pairs | WHO classes |
|---|---|---:|---:|---:|---:|
| 02a | Harmonized lineage groups | 803,628 | 2,832 | 23,196 | — |
| 02b | Retained groups mapped to a historical WHO class | 697,116 | 1,750 | 19,129 | 13 |
| 02c | Retained country–variant pairs whose peak daily share is `> 1%` | 564,232 | 1,613 | 12,078 | 13 |
| 03 | Retained country–variant pairs with at least 7 consecutive calendar days at `share > 1%`, starting at the first qualifying window | 206,994 | 566 | 1,753 | 10 |

The 02c output is `code/variants_who.rds`; the 03 output is
`code/variants_who2.rds`. Step 03 removes all observations before the first day
of the first qualifying seven-day window, then keeps the remaining time series.
Later observations may have `share <= 1%`. A missing date interrupts a
consecutive run.

After step 03, Lambda, Theta, and Zeta are absent because none of their
country–variant pairs pass the seven-day rule.

The 03 output contains only these seven core columns:
`country`, `date`, `variant`, `who_variant`, `denominator`, `numerator`, and
`share`.

## WHO classes and approximate dates

The workflow defines 13 historical WHO classes. For classes retained by step
03, the date range below runs from the earliest first qualifying window to the
latest qualifying-window end across countries. For excluded classes, it shows
their full observed range in the 02 data. These are approximate dataset ranges,
not official WHO designation dates.

| WHO class | Pango root(s) | Approximate dates in these data | Present after 03? |
|---|---|---|---|
| Alpha | B.1.1.7 | 2020-11-08 to 2021-09-03 | Yes |
| Beta | B.1.351 | 2021-01-18 to 2021-07-28 | Yes |
| Gamma | P.1 / B.1.1.28.1 | 2021-02-05 to 2021-07-29 | Yes |
| Delta | B.1.617.2 | 2021-04-12 to 2022-01-28 | Yes |
| Epsilon | B.1.427 and B.1.429 | 2020-11-28 to 2021-05-19 | Yes |
| Zeta | P.2 / B.1.1.28.2 | Observed 2020-10-11 to 2021-11-29; no qualifying window | No |
| Eta | B.1.525 | 2021-02-15 to 2021-05-08 | Yes |
| Theta | P.3 / B.1.1.28.3 | Observed 2021-02-25 to 2021-09-01; no qualifying window | No |
| Iota | B.1.526 | 2021-01-25 to 2021-06-29 | Yes |
| Kappa | B.1.617.1 | 2021-03-22 to 2021-04-10 | Yes |
| Lambda | C.37 / B.1.1.1.37 | Observed 2021-02-01 to 2021-10-12; no qualifying window | No |
| Mu | B.1.621 | 2021-04-25 to 2021-07-23 | Yes |
| Omicron | B.1.1.529 and descendants | 2021-12-01 to 2024-09-09 | Yes |

## Counting definitions

- **Unique Pango variants** counts distinct values of `variant` across all countries.
- **Country–variant pairs** counts distinct `country + variant` combinations.
- **WHO classes** counts distinct non-missing values of `who_variant`.

## Re-run note

The current `02_group_variants.R` source also contains a seven-day filter, while
the saved `variants_who.rds` reflects the earlier 02c output shown above. Running
the current 02 script as-is may therefore overwrite `variants_who.rds` with an
already seven-day-filtered dataset, making step 03 redundant. For the documented
two-script workflow, apply the seven-day filter only in `03_filter_7day.R`.
