# Old GLMS David wrote
####################
## BIRD OBS MODELS #
####################
bird.obs.nb <- glmer.nb(count~treatment+ (1|site), 
                 	data=bird.obs) 
sim.bird.obs.nb <- simulateResiduals(fittedModel = bird.obs.nb, n = 250)
plot(sim.bird.obs.nb)

bird.obs.nb.zi <- glmmTMB(count~treatment+ (1|site), 
                      family=nbinom1, 
                      zi=~1, 
                      data=bird.obs)
sim.bird.obs.nb.zi <- simulateResiduals(fittedModel = bird.obs.nb.zi, n = 250)
plot(sim.bird.obs.nb.zi)
Anova(bird.obs.nb.zi)

bird.obs.pois.zi <- glmmTMB(count~treatment+ (1|site), 
                      family=poisson, 
                      zi=~1, 
                      data=bird.obs)
sim.bird.obs.pois.zi <- simulateResiduals(fittedModel = bird.obs.pois.zi, n = 250)
plot(sim.bird.obs.pois.zi)
Anova(bird.obs.pois.zi)
####################
# BIRD RICH MODELS #
####################
# GLMM
glmm.bird.rich <- glmer(richness ~ treatment + (1|sites),
												data=bird.rich, family=poisson)
sim.glmm.bird.rich <- simulateResiduals(fittedModel = glmm.bird.rich, n = 250)
plot(sim.glmm.bird.rich) # fails dispersion

# OD GLMM
bird.rich$obs <- 1:length(bird.rich$sites)
glmm.ovd.bird.rich <- glmer(richness ~ treatment + (1|sites) + (1|obs),
												data=bird.rich, family=poisson)
sim.glmm.ovd.bird.rich <- simulateResiduals(fittedModel = glmm.ovd.bird.rich, n = 250)
plot(sim.glmm.ovd.bird.rich)
hist(sim.glmm.ovd.bird.rich)

Anova(glmm.ovd.bird.rich)
####################
# SEED OBS MODELS #
####################
# Poisson zero inflated
poisZI.seeds.obs <- glmmTMB(SEEDS~TREATMENT+
                        (1|BLOCK), 
                      family=poisson, 
                      zi=~1, 
                      data=seed.obs)
sim.poisZI.seeds.obs <- simulateResiduals(fittedModel = poisZI.seeds.obs, n = 250)
plot(sim.poisZI.seeds.obs)
Anova(poisZI.seeds.obs)

####################
# SEED RICH MODELS #
####################

poisZI.seed.rich <- glmmTMB(richness~treatment+
                        (1|sites), 
                      family=poisson, 
                      zi=~1, 
                      data=seed.rich)
sim.poisZI.seed.rich <- simulateResiduals(fittedModel = poisZI.seed.rich, n = 250)
plot(sim.poisZI.seed.rich)
Anova(poisZI.seed.rich)
