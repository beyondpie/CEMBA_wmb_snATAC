#!/bin/bash

# * mm10 lift to hg38
liftOver whole.mouse.brain.cCREs.bed mm10ToHg38.over.chain.gz -minMatch=0.5 peak.conserve0.5.bed peak.unMapped0.5

# * hg38 lift back to mm10
liftOver peak.conserve0.5.bed hg38ToMm10.over.chain.gz -minMatch=0.5 peak.conserve0.5.reciprocal0.5.bed peak.conserver0.5.unMapped0.5 

