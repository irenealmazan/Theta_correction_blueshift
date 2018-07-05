function [estimates, model] = fitpc(Im, Icoh, lx, ly)

start_point = [lx ly];

model = @gaussblur;

estimates = fminsearch(model, start_point);

    function [sse, FittedCurve] = gaussblur(params)

        sigX = params(1);
        sigY = params(2);    
        [sz1 sz2 sz3] = size(Im);
        [X Y]= meshgrid([-sz1/2:sz1/2-1], [-sz2/2:sz2/2-1]);
        
        g= exp( -(X.^2) / (2*sigX^2) -(Y.^2) / (2*sigY^2));
        
        for ii=1:sz3
            FittedCurve(:,:,ii) = fftshift(ifftn( fftn(fftshift(Icoh(:,:,ii))).*fftn(fftshift(g)) ));
        end
            
        ErrorVector = FittedCurve - Im;
        
        
        %Im=Im/mean(Im(:));
        %Icoh = Icoh/mean(Icoh(:));
        %g=g/mean(g(:));
        %g=g*10000;
        
        %fIm = (fftn(Im));
        %fIcoh = (fftn(Icoh));
        %fg = (fftn(fftshift(g)));
        %fIm(1)=0; fIcoh(1)=0; fg(1)=0;
        %fIm=fftshift(fIm); fIcoh=fftshift(fIcoh); fg=fftshift(fg);
        
        %FittedCurve = fftn(fftshift(Icoh)).*fftn(fftshift(g));
        %ErrorVector = FittedCurve - ifftn(Im);
         
        sse = sum(ErrorVector(:).^2);

    end

end