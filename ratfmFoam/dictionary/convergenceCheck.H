// Check convergence
if (
	maxResidualp < convergenceCriterion
	&& maxResidualk1 < convergenceCriterion
	&& maxResidualk2 < convergenceCriterion
	&& maxResiduale1 < convergenceCriterion
	&& maxResiduale2 < convergenceCriterion
	&& maxResidualt < convergenceCriterion
	)
{
    Info<< "reached convergence criterion: " << convergenceCriterion << endl;
    runTime.writeAndEnd();
    Info<< "latestTime = " << runTime.timeName() << endl;
}

