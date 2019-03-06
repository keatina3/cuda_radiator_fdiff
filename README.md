<article data-history-node-id="1724" role="article" about="/node/1724" class="node node--type-assignment node--view-mode-full clearfix">
  <header>
    
        
<div class="node__meta">
<article typeof="schema:Person" about="/user/18" class="profile">
  <a href="/blog/18">View recent blog entries</a></article>

<span>
Submitted by <span class="field field--name-uid field--type-entity-reference field--label-hidden"><a title="View user profile." href="/user/18" lang="" about="/user/18" typeof="schema:Person" property="schema:name" datatype="" class="username">jose</a></span>
 on <span class="field field--name-created field--type-created field--label-hidden">Fri, 2019/03/01 - 10:06</span>
</span>
        
</div>
</header>
  <div class="node__content clearfix">
    
  <div class="field field--name-field-course-reference field--type-entity-reference field--label-inline">
<div class="field__label">Course Reference</div>
<div class="field__item"><a href="/node/5" hreflang="en">5615</a></div>
</div>

  <div class="field field--name-field-assignment-number field--type-integer field--label-inline">
<div class="field__label">Assignment Number</div>
<div class="field__item">2</div>
</div>

  <div class="field field--name-field-due-date field--type-datetime field--label-inline">
<div class="field__label">Due Date</div>
<div class="field__item"><time datetime="2019-03-20T17:00:00Z" class="datetime">Wed, 2019/03/20 - 17:00</time>
</div>
</div>

<div class="clearfix text-formatted field field--name-body field--type-text-with-summary field--label-hidden field__item"><h2>5615-02 - Cylindrical Radiator Finite Differences model</h2>

<p>The goal of this assignment is to model the propagation of heat inside a cylindrical radiator.</p>

<p>===================================================================================================================</p>

<p><strong>Task 1 - calculation in the cpu</strong></p>

<p>Write C or C++ code in .h/.c or .hpp/.cpp files that allocates two <strong>floating (not double)</strong> point matrices of size <em>n</em> x <em>m</em> (both with a default value of 32 and that can be specified by passing command line arguments -n and -m, respectively). The default number of iterations will be 10, but will have to be specified with a -p number command line argument option.</p>

<p>The boundary conditions will be:</p>

<ul><li>in column 0, matrix[i][0] = 1.00*(float)(i+1)/(float)(n) (so the values range between 1.00/(float)(n) and 1.00)</li>
	<li>in column 1, matrix[i][1] = 0.80*(float)(i+1)/(float)(n) (so the values range between 0.80/(float)(n) and 0.80)</li>
</ul><p>Those values will remain constant and thus columns 0 and 1 do not need to be calculated on each time step. Initial conditions will be 0 in any other position of the matrices that isn't in columns 0 or 1.</p>

<p>Propagation of heat happens only in rows (so there is no propagation vertically) and is directional (in the sense that water in the radiator flows mostly towards "the right"), so you must apply the following weights:</p>

<p>(*nextMatrix)[ui][uj]= ( (1.9*(*previousMatrix)[ui][uj-2])+(1.5*(*previousMatrix)[ui][uj-1])+ (*previousMatrix)[ui][uj ]+ (0.5*(*previousMatrix)[ui][uj+1])+(0.1*(*previousMatrix)[ui][uj+2]) ); (*nextMatrix)[ui][uj]/=(float)(5.0);</p>

<p>The radiator is horizontally cylindrical (and a bit of heat propagates leftwise as well) and each row (pipe) is a cycle, so, for example, to compute the positions:</p>

<ul><li>new [ui][m-2], you will need the values of old [ui][m-4],[ui][m-3],[ui][m-2],[ui][m-1] and [ui][0]</li>
	<li>new [ui][m-1], you will need the values of old [ui][m-3],[ui][m-2],[ui][m-1],[ui][0] and [ui][1]</li>
</ul><p>Add an command line option (say, "-a"), so that, after the designated number of heat propagation time steps has been concluded, the average temperature for each row of the radiator gets calculated (this represents the thermostat that would be used to stop heating once it has reached a certain level).</p>

<p><strong>Task 2 - parallel implementation</strong></p>

<p>From the code in the host, add cuda code in .h and .cu files that implements the same model in the gpu. To make things easier, you can make the following assumptions:</p>

<p>* The <em>final</em> <em>n</em> and <em>m</em> sizes are expected to be a multiple of 32 (rectangular cases, in which m != n, must still work, though) - it is recommended to start with smaller problem and block sizes (say, 5), for easier debugging.</p>

<p>* Feel free to use atomics, as well as shared, texture or surface memory, as well as any other CUDA feature (but not other libraries other than CUDA) to speed up the calculation. Feel free to check and reuse the provided source code (for example cudaEvents, atomic, sharedMemory, etc).</p>

<p>* You can organize the grid in whatever way you want.</p>

<p>You will need to copy back the results to the RAM; compare them with the values obtained in the CPU and report any mismatches larger than 1.E-5. Add code to track the amount of time (use cuda events for the gpu parts) that each one of the steps takes (compute on the cpu, allocation on the gpu, transfer to the gpu, compute on the gpu, calculation of the averages, transfer back to the ram) takes, and add a command line argument (-t) to display both the CPU and GPU timings next to each other, as well as the other steps.</p>

<p><strong>Task 3 - performance improvement</strong></p>

<p>Test the code for different sizes, use <em>n=15360</em>,<em>m=15360</em> and <em>p=1000</em> as a reference. Try as well different numbers of threads per block, and different numbers in each one of the x and y directions, if you are using a 2d grid. Calculate the speedups (CPU vs GPU) and precisions compared to the CPU versions, first for compute, and second, including as well the memory allocations and transfers.</p>

<p>Comment on how the new reduce that you have implemented to calculate the average temperature per row performs compared to the reduce operation that you implemented for the first assignment.</p>

<p><strong>Note:</strong> In this assignment, for decent CPU code, you can expect reasonable (&gt;6) but not particulaly large (&lt;60) speedups (for single precision, at least!)- however, the lower that the cuda execution times are, the more marks that the assignment will receive.</p>

<p><strong>Task 4 - double precision version</strong> Port the code to a double precision format and compare the times, speedups and mismatches, if any.</p>

<p>===================================================================================================================</p>

<p>Submit a tar ball with your source code files (including a working Makefile), speedup graphs and a writeup of what you did and any observations you have made on the behaviour and performance of your code, as well as problems that you came across while writing the assignment.</p>

<p><strong>Note:</strong> Do not forget to set up your PATH and LD_LIBRARY_PATH variables in your .bashrc, otherwise nvcc won't work! You have the instructions in the slide 65 of the pdf of the course.</p>

<p><strong>Note:</strong> Marks will be deducted for tarbombing. http://en.wikipedia.org/wiki/Tar_%28computing%29#Tarbomb</p>

<p><strong>Note:</strong> Extra marks will be given for separating the C/C++ and cuda code. You can find examples on how to do that on the "makefileCpp" and "makefileExternC" sample code.</p>

<p><strong>Note:</strong> Remember that the code must work for non-square systems too, even when, for this assignment, we are using square ones for benchmarking. You can run use <strong>cuda-memcheck</strong> to test that, as in: cuda-memcheck ./my_exec</p>

<p><strong>Note:</strong> When you are benchmarking the performance of your code, you can check the current load on cuda01 with the <em>nvidia-smi</em> command. ===========================================================================================================================</p></div>
      
  </div>
</article>

