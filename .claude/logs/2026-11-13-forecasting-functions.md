# 2026-11-13: forecasting functions (v0.4.0)

For MPMS (code/model-julia), whose next-period consumption stage needs this period's portfolio.

- `ModelSolution` gains type parameter `FF` and field `Ffct` (default `nothing`).
- Hook `forecasts_from_policy(state, policy, mdl)` (stubs.jl), default `(;)`. Non-empty: `main_setup`
  builds `Ffct` (warm start from `initobj.Ffct`); `forecast_next(trans, obj, forecasts)` merges the
  values into every `next` point.
- STPFIMethod, fixed-point steady state, simulation (`calc_eqm_from_state`): forecasts from the
  policy being evaluated. TFIMethod: stored `Ffct` (frozen during the solve), updated and damped
  like `Vfct`; its distance (`distF`) counts toward convergence and is logged (F line, CSV columns).
  TFIMethod needs equations that depend on the policy outside the forecasts (docstring).
- `FunctionValuesAtPoint` gains `ffct` (4-arg constructor keeps `nothing`).
- Tests: test/test_forecasts.jl (NGM with the Euler kernel using current C as a forecast; TFI and
  STPFI reproduce NGM; FP steady state; zero allocations). 824/824 pass.
- Not registered yet (user approval pending).
