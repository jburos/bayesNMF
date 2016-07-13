## from : https://gist.github.com/danstowell/151a3b68ceb9c5f4a51f
data {
	int<lower=1> N;  // num frames
	int<lower=1> M;  // num bins
	int<lower=1> K;  // num factors
	real<lower=0> Wconc; // concentration parameter on template priors - a large number means relative certainty, stick close to the prior
	real<lower=0> Hconc; // concentration parameter on activations - a large number encourages sparsity
	matrix<lower=0>[M,K] Winit; // prior spectral templates -- these should each be different (or else it'd be unidentifiable). avoid zeroes too.
	matrix<lower=0>[M,N] X; // spectrogram
}
transformed data {
	vector<lower=0>[M*N] Xnorm;
	vector<lower=0>[M] Winitnorm[K];

	Xnorm = to_vector(X / sum(X)); // spectrogram, normalised and flattened
	for (k in 1:K) {
		Winitnorm[k] = col(Winit, k) * Wconc / sum(col(Winit, k)); // Winit, normalised and pre-multiplied by concn
	}
}
parameters {
	simplex[M] W[K]; // inferred spectral templates
	matrix<lower=0>[K,N] H; // activations
	real<lower=0> Xconc; // overall concentration of spectrogram, which means how closely it sticks to the dotproduct, i.e. relates to noise level
}
transformed parameters {
	vector[M*N] Xest;
	matrix[K,M] Wmat; // W, but indexed backwards and transposed later, since it seems that's the only way to assign the simplices into the matrix columns

	for (k in 1:K) {
		Wmat[k] = W[k]';
	}
	Xest = to_vector(Wmat' * H);
}
model {
	for (k in 1:K) {
		W[k] ~ dirichlet(Winitnorm[k]); // our supplied spectral templates are used as centres for the priors on the spectral templates. each tpl has its own dirichlet.
	}
	to_vector(H) ~ exponential(Hconc); // a prior on the activations

	Xconc ~ lognormal(100, 10); // a high value of Xconc means low SNR (a more semantic parametrisation would be nice...)
	Xnorm ~ dirichlet(Xest * Xconc); // our input data is modelled as a multinomial sampled from this overall mega-dirichlet
}
