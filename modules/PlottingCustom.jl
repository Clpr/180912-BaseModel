"""
    PlottingCustom

Customized methods of plotting, using PyPlot
"""
module PlottingCustom
    using PyPlot
    export QuickSubplot, QuickSubplot_Custom1
# ==============================================================================
    function QuickSubplot(ROW::Int64,COL::Int64,POS::Int64,
        X,Y,Y_twin,
        TITLE::String,XLAB::String,YLAB::String,LEGEND_twin::Array{String,1},
        XLIM::Array,YLIM::Array,YLIM_twin::Array,YLAB_twin::String  )
        # --------
        # quick-style macros for m*n subplot panels
        obj = subplot(ROW,COL,POS)
                plot(X,Y);title(TITLE);xlabel(XLAB);ylabel(YLAB);grid(true);
                ylim(YLIM)
                twinx();plot(X,Y_twin,"--r");ylabel(YLAB_twin);
                legend(LEGEND_twin,loc="best")
                xlim(XLIM)
                ylim(YLIM_twin)
        return obj
    end
    function QuickSubplot_Custom1(ROW::Int64,COL::Int64,POS::Int64,X,Y,Y_twin,TITLE,XLAB,LEGEND_twin,IDXRANGE,TRANGE; DELTA=0.1)
        tmp_Yupper = findmax(Y[IDXRANGE[1]:IDXRANGE[2]])[1]
        tmp_Ylower = findmin(Y[IDXRANGE[1]:IDXRANGE[2]])[1]
        tmp_Ydelta = (tmp_Yupper - tmp_Ylower) * DELTA
        tmp_Yupper_twin = findmax(Y_twin[IDXRANGE[1]:IDXRANGE[2]])[1]
        tmp_Ylower_twin = findmin(Y_twin[IDXRANGE[1]:IDXRANGE[2]])[1]
        tmp_Ydelta_twin = (tmp_Yupper_twin - tmp_Ylower_twin) * DELTA

        obj = QuickSubplot(ROW,COL,POS,  X,Y,Y_twin,
                TITLE,XLAB,TITLE,LEGEND_twin,
                TRANGE  ,
                [tmp_Ylower - tmp_Ydelta  , tmp_Yupper + tmp_Ydelta  ],
                [tmp_Ylower_twin - tmp_Ydelta_twin  , tmp_Yupper_twin + tmp_Ydelta_twin  ],
                LEGEND_twin[1]   )
        return obj
    end





























# ==============================================================================
end  # module ends
#
