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

## WHO classes and variant start dates

Step 03 retains the 10 WHO classes shown below. For each unique Pango variant,
its start date is the earliest date across countries on which its retained time
series begins—that is, the first day of its first qualifying seven-day window.
These are dataset-derived dates, not official WHO designation dates.

On 5 May 2023, WHO declared that COVID-19 no longer constituted a public health
emergency of international concern ([PHEIC](https://www.who.int/europe/news/item/05-05-2023-statement-on-the-fifteenth-meeting-of-the-international-health-regulations-%282005%29-emergency-committee-regarding-the-coronavirus-disease-%28covid-19%29-pandemic)).
In this dataset, 153 variants—all Omicron—started after this date. No variant
started exactly on 5 May 2023.

| WHO class | Variant count | Earliest start | Latest start |
|---|---:|---|---|
| Alpha | 1 | 2020-11-08 | 2020-11-08 |
| Epsilon | 2 | 2020-11-28 | 2020-12-09 |
| Beta | 3 | 2021-01-18 | 2021-06-19 |
| Iota | 1 | 2021-01-25 | 2021-01-25 |
| Gamma | 6 | 2021-02-05 | 2021-05-24 |
| Eta | 1 | 2021-02-15 | 2021-02-15 |
| Kappa | 1 | 2021-03-22 | 2021-03-22 |
| Delta | 95 | 2021-04-12 | 2021-11-29 |
| Mu | 1 | 2021-04-25 | 2021-04-25 |
| Omicron (before 2023-05-05) | 302 | 2021-12-01 | 2023-05-03 |
| $\color{red}{\textsf{Omicron\ (after\ 2023-05-05)}}$ | $\color{red}{153}$ | $\color{red}{\textsf{2023-05-08}}$ | $\color{red}{\textsf{2024-08-11}}$ |

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
