# Krüger 5 nm authority attempt `be95cd4`: numerical refusal

This attempt is not an experimental or calibration result. It stopped after 23 accepted steps at
1.963905945869131 s physical time because one slightly downward diffuse-neutral ray required more
than the then-authorized 1,024 exact periodic-cell traversals.

The saved checkpoint reproduced the ray exactly:

- origin: `[0.108438417364833, 0.004037909820028365, 2.134689932283348]`
- direction: `[-0.00391740694902086, -0.9999561259207863, -0.008508828138368412]`
- exact outcome: hard hit on face 805 after 1,110 wraps
- exact hit position: `[0.021500688484835793, 0.012339212119627343, 1.94585629732885]`
- conservative lower-domain integrity horizon: 12,556 wraps

Thus the ray is finite and physical; the fixed/upward-only replay horizon was incomplete. Commit
`a15634f` generalizes the authority to a geometry-derived vertical-domain horizon. Upward rays may
hard-hit or escape through the open top. Downward rays must hard-hit before the derived lower-domain
horizon. Exactly horizontal rays and downward rays without a hard hit still refuse.

Artifact checksums:

- `failed_audit.json`: `24f3463c1e79e5dda5fe8de7fd673c651a701d70c283d24e4eb26ba37ed11838`
- `failed_run.log`: `d95b428db12ffa26a907cb8e98e3ae707c141a573bd43d0f67af1ce72b2c45b6`
- `failed_checkpoint.npz`: `6ae5daee01e2211222bc28d6901ae830d8e19a1f0081b99a69feecfc9dcdb53a`

The authoritative retry restarts from zero; this checkpoint is preserved only as regression
evidence and must never be promoted or resumed into an authority result.
