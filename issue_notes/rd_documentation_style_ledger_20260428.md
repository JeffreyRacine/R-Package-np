# Rd Documentation Style Ledger

Date: 2026-04-28

Status: provisional until the first live pilot family passes validation.

## Purpose

This ledger records reusable `.Rd` documentation patterns for the release
documentation campaign. It is a style contract, not a status tracker.

## Provisional Argument-Group Names

- `Data And Bandwidth Inputs`
- `Formula Interface`
- `Bandwidth Selection Controls`
- `Continuous Scale-Factor Search Initialization`
- `Discrete And Categorical Search Initialization`
- `Kernel Type Controls`
- `Continuous Kernel Support Controls`
- `Local-Polynomial And Degree Search Controls`
- `Least-Squares Quadrature Controls`
- `Numerical Search And Tolerance Controls`
- `Plot And Interval Controls`
- `Returned Values`
- `Compatibility And Additional Arguments`

Use these names in details-level argument guides or as planning labels. Do not
put `\subsection{}` inside `\arguments{}`.

## Provisional Reader QA Rubric

For touched public pages, rendered help should answer these questions quickly:

1. What is the main call pattern?
2. Should the user call the estimator or the `*bw` selector?
3. Where are bandwidth/search controls documented?
4. Which controls define the estimator and which define bandwidth selection?
5. What does `...` pass through?
6. Where should the user go next: estimator, `*bw`, kernel/options page,
   vignette, or Gallery?

Rendered output is the review surface. Source neatness is not enough.

## Provisional Estimator-To-Bandwidth Wording

Use this pattern only after adapting it to the function family:

> If `bws` is omitted, `FUNCTION()` computes bandwidths by calling
> `FUNCTIONbw()`. Bandwidth-selection controls supplied through `...` are passed
> to `FUNCTIONbw()`; see `?FUNCTIONbw` for the full set of bandwidth, kernel,
> support, search, local-polynomial, quadrature, and scale-factor controls.

## Provisional `bws` Wording

> An optional bandwidth object, usually returned by `FUNCTIONbw()`. If supplied,
> it is used directly. If omitted, `FUNCTION()` computes bandwidths internally by
> calling `FUNCTIONbw()`.

## Provisional `...` Wording

> Additional arguments. When bandwidths are not supplied, bandwidth-selection
> controls in `...` are passed to `FUNCTIONbw()`. Other arguments are used by the
> estimator or related methods as documented for this page.

## Provisional Gallery Wording

To be filled only after confirming the canonical Gallery URL/text from existing
vignettes. Do not invent a new URL or wording during `.Rd` edits.

## Do Not Use

- Do not put `\subsection{}` blocks inside `\arguments{}`; touched-page
  validation showed that this can make following `\item{}` entries parse
  incorrectly.
- Do not create fake group headings with `\item{Group Name}{...}` inside
  `\arguments{}`.
- Do not reorder R function formals during the first-pass documentation
  campaign.
- Do not reorder `.Rd` `\usage{}` for visual tidiness.
- Do not copy the full `*bw` argument universe into estimator pages.
- Do not over-section moderate pages; too many headings can make help harder to
  scan.
