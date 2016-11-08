classdef test_phase_surrogates < matlab.unittest.TestCase

	% Some unit tests to verify properties of the phase randomization process

	methods(Test)
		function test_odd_samples(self)
			cmat = [1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1]; % Noise covariance
			x = mvnrnd([0 1 2],cmat,1001);

			y = ROInets.generate_phase_surrogates(x);
			
			% Check means are preserved
			self.verifyEqual(mean(x),mean(y),'AbsTol',1e-5);

			% Check STD is preserved
			self.verifyEqual(std(x,1),std(y,1),'AbsTol',1e-5);

			% Check correlations across channels are destroyed
			x_corr = triu(cov(y),1);
			self.verifyEqual(x_corr(:),zeros(numel(x_corr),1),'AbsTol',0.15); % They won't be totally destroyed, but should be nowhere near 0.5

			% Check power spectra are the same
			self.verifyEqual(abs(fft(x)).^2,abs(fft(y)).^2,'RelTol',1e-5); %

		end

		function test_even_samples(self)
			cmat = [1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1]; % Noise covariance
			x = mvnrnd([3 1 2],cmat,1000);

			y = ROInets.generate_phase_surrogates(x);
			
			% Check means are preserved
			self.verifyEqual(mean(x),mean(y),'AbsTol',1e-5);

			% Check STD is preserved
			self.verifyEqual(std(x,1),std(y,1),'AbsTol',1e-5);

			% Check correlations across channels are destroyed
			x_corr = triu(cov(y),1);
			self.verifyEqual(x_corr(:),zeros(numel(x_corr),1),'AbsTol',0.15); % They won't be totally destroyed, but should be nowhere near 0.5

			% Check power spectra are the same
			self.verifyEqual(abs(fft(x)).^2,abs(fft(y)).^2,'RelTol',1e-5); 

		end

		function test_preserve_correlation(self)
			cmat = [1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1]; % Noise covariance
			x = mvnrnd([1 1 2],cmat,1000);

			y = ROInets.generate_phase_surrogates(x,true);
			
			% Check means are preserved
			self.verifyEqual(mean(x),mean(y),'AbsTol',1e-5);

			% Check STD is preserved
			self.verifyEqual(std(x,1),std(y,1),'AbsTol',1e-5);

			% Check correlations across channels are PRESERVED now
			self.verifyEqual(cov(x),cov(y),'AbsTol',1e-5); 

			% Check power spectra are the same
			self.verifyEqual(abs(fft(x)).^2,abs(fft(y)).^2,'RelTol',1e-5); 

		end


		
	end
end

