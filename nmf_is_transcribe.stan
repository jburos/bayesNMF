## from : https://gist.github.com/danstowell/151a3b68ceb9c5f4a51f

data {
	int<lower=1> N;  // num patients
	int<lower=1> M;  // num trinucleotides
	int<lower=1> K;  // num signatures
	real<lower=0> Wconc;        // concentration parameter on template priors - a large number means relative certainty, stick close to the prior
	real<lower=0> Hconc;        // concentration parameter on activations - a large number encourages sparsity
	matrix<lower=0>[M,K] Winit; // prior spectral templates -- these should each be different (or else it'd be unidentifiable). avoid zeroes too.
	matrix<lower=0>[M,N] X;     // spectrogram
}
transformed data {
	vector<lower=0>[M*N] Xnorm;
	vector<lower=0>[M] Winitnorm[K];

	Xnorm = to_vector(X); // spectrogram, flattened, but in this case NOT normalised (since not PLCA)
	for (k in 1:K) {
		Winitnorm[k] = col(Winit, k) * Wconc / sum(col(Winit, k)); // Winit, normalised and pre-multiplied by concn
		//for (m in 1:M) {
		//  print("k = ", k, "; m = ", m, "; val = ", Winitnorm[k,m]);
		//}
	}
}
parameters {
	simplex[M] W[K];        // inferred spectral templates
	matrix<lower=0>[K,N] H; // activations
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
    // our supplied spectral templates are used as centres for the priors on the spectral templates. each tpl has its own dirichlet.
		W[k] ~ dirichlet(Winitnorm[k]); 
	}
	to_vector(H) ~ exponential(Hconc); // a prior on the activations
	
  // IS-NMF -- the factor of 1/2 here comes from X being assumed a power-spectrogram where the underlying random variables are complex unit-norm
	Xnorm ~ normal(0, 0.5 * Xest); 
}
generated quantities {
  matrix[M,N] Xhat;
  vector[M*N] Xhat_pre;
  
  for (i in 1:(M*N)) {
    Xhat_pre[i] = normal_rng(0, 0.5 * Xest[i]);
  }
  for (m in 1:M) { ## rows
    for (n in 1:N) { ## columns
      Xhat[m,n] = Xhat_pre[(n-1)*M+m];
    }
  }
}
