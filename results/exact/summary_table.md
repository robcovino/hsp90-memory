| model | description | k | max log L | BIC | log Z | Δlog Z | held-out log L (CV) | note |
|---|---|---|---|---|---|---|---|---|
| modal3 | chain of 3 modes | 12 | -98360.5 | 196842.7 | -98416.3 | -407.8 | — |  |
| modal4 | chain of 4 modes | 17 | -98133.7 | 196439.9 | -98206.3 | -197.8 | -98177.2 (cv_v4) |  |
| modal5 | chain of 5 modes | 22 | -98092.3 | 196407.9 | -98181.8 | -173.4 | -98146.6 (cv_v5) |  |
| modal6 | chain of 6 modes | 27 | -98058.0 | 196390.1 | -98162.5 | -154.1 | -98115.8 (cv_v4) |  |
| modal7 | chain of 7 modes | 32 | -98013.2 | 196351.2 | -98157.2 | -148.7 | -98087.9 (cv_v5) | multimodal posterior: log Z of the dominant sampled mode |
| modal3_ou | chain 3 + OU on kOC | 14 | -98217.0 | 196576.1 | — | — | — |  |
| modal4_ou | chain 4 + OU on kOC (non-eq.) | 19 | -98056.6 | 196306.0 | -98134.0 | -125.6 | -98101.2 (cv_v4) |  |
| modal4_oueq | chain 4 + OU on kOC (eq.) | 19 | -98056.9 | 196306.6 | -98134.3 | -125.8 | -98100.3 (cv_v5) |  |
| modal4_ou2 | chain 4 + OU on kOC and kCO | 20 | -98056.4 | 196315.7 | — | — | — |  |
| modal5_oueqF | chain 5 + OU (eq.) | 24 | -97998.6 | 196240.9 | -98127.8 | -119.4 | -98052.5 (cv_v6) | log Z from an inferior peak (ML -98034.7 at the time) |
| ring4 | ring of 4 modes | 19 | -98050.0 | 196292.8 | -98125.6 | -117.1 | -98097.3 (cv_v6) |  |
| ring4_noDB | ring of 4, cycle not balanced | 20 | -98050.2 | 196303.4 | — | — | — |  |
| ringL1 | ring 4 + 1 discrete kOC ladder mode | 24 | -98007.7 | 196259.0 | -98100.6 | -92.1 | -98064.4 (cv_v5) |  |
| ringL2 | ring 4 + 2-step discrete kOC ladder | 29 | -97977.3 | 196249.0 | -98084.9 | -76.4 | -98037.9 (cv_v5) |  |
| ringS2 | ring 4 + 2 discrete kOC branches | 29 | -97976.3 | 196246.9 | -98084.3 | -75.8 | -98042.0 (cv_v6) |  |
| ringF1 | ring 4 + second frozen state | 24 | -97996.3 | 196236.2 | -98080.3 | -71.8 | -98040.0 (cv_v6) |  |
| ring4_oueq | ring 4 + OU on kOC (eq.) | 21 | -97971.2 | 196155.4 | -98052.1 | -43.7 | -98017.1 (cv_v6) |  |
| ring4_noDB_oueq | ring 4 + OU, cycle not balanced | 22 | -97971.1 | 196165.5 | — | — | -98022.8 (cv_v5) |  |
| ringF1_oueq | ring 4 + second frozen state + OU (eq.) | 26 | -97918.9 | 196101.6 | -98008.4 | 0.0 | -97964.1 (cv_v6) |  |
| full5_oueq | all 10 connections of 5 modes + OU (eq.) | 36 | -97909.2 | 196183.7 | — | — | -97957.0 (cv_v6) | likelihood ceiling for 5 modes |
