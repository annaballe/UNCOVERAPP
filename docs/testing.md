# Test e validazione

Unit test (testthat):
- tests/testthat/ e inst/tests/testthat/ contengono test per funzioni di libreria.
- Esecuzione:
```r
devtools::test()
```

Test dell’app Shiny:
- inst/shiny-dir/tests/ contiene script e output attesi (JSON/PNG) per scenari base.
- Suggerito l’uso di shinytest2 per test end-to-end (opzionale).
