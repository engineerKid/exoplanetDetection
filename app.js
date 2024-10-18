// Define the openTab function
function openTab(evt, tabName) {
    let i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
}

// Ensure the first tab is open by default when the page loads
document.addEventListener('DOMContentLoaded', function () {
    // Automatically click the first tab to display it
    const firstTab = document.querySelector('.tablinks');
    if (firstTab) {
        firstTab.click();
    }
});

// Info box - Toggle Collapse Functionality
function toggleCollapse(id) {
    const content = document.getElementById(id);
    content.style.display = content.style.display === 'none' ? 'block' : 'none';
}

// Ensure the collapsible boxes are initially hidden
document.querySelectorAll('.content-box').forEach(box => box.style.display = 'none');

// Existing DOMContentLoaded event listener
document.addEventListener('DOMContentLoaded', function () {

    // Use the dynamically injected RA_CENTER, DEC_CENTER, and SEARCH_RADIUS values
    const raMin = RA_CENTER - SEARCH_RADIUS;
    const raMax = RA_CENTER + SEARCH_RADIUS;
    const decMin = DEC_CENTER - SEARCH_RADIUS;
    const decMax = DEC_CENTER + SEARCH_RADIUS;

    fetch('/get_star_data')
        .then(response => response.json())
        .then(data => {
            const starData = data;

            // Creates two plots: (1)sky region plot and (2)H-R diagram

            /**
             * ********************************
             * #1: Code for the sky region plot
             * ********************************* 
             */

            // Data for sky region plot
            const ra = starData.map(star => star.ra);
            const dec = starData.map(star => star.dec);
            const tmag = starData.map(star => star.Tmag);
            const hasExoplanet = starData.map(star => star.c_has_exoplanet ? 'Exoplanet' : 'Star');

            // Prepare the sky region scatter plot
            const skyTrace = {
                x: ra,  // RA
                y: dec,  // DEC
                mode: 'markers',
                marker: {
                    size: tmag.map(t => 16 - t),  // Size inversely proportional to Tmag
                    color: hasExoplanet.map(e => e === 'Exoplanet' ? 'red' : 'blue'), // exoplanet host star denoted be red
                    opacity: 0.8,  // Higher opacity for better visibility
                },
                text: hasExoplanet,
                hoverinfo: 'text'
            };

            // To restrict panning/zooming within the set ra/dec boundary
            const skyLayout = {
                title: 'Sky Region: Stars and Exoplanets',
                xaxis: {
                    title: 'Right Ascension (RA)',
                    range: [raMin, raMax],  // Dynamic RA boundary
                    autorange: false,
                    constrain: 'domain',
                },
                yaxis: {
                    title: 'Declination (DEC)',
                    range: [decMin, decMax],  // Dynamic Dec boundary
                    autorange: false,
                    constrain: 'domain',
                },
                showlegend: false,
                scrollZoom: true,  // Allow zooming within the boundary
                dragmode: 'zoom',  // Allow panning within the boundary
            };

            // Plot the sky region  
            Plotly.newPlot('sky-plot', [skyTrace], skyLayout);
           
            /**
             * ******************************
             * #2: Code for H-R diagram plot 
             * ******************************
             */
            // Data for H-R diagram
            const teff = starData.map(star => star.Teff);
            const lum = starData.map(star => star.lum);

             // Prepare the H-R diagram scatter plot
             const hrTrace = {
                x: teff,
                y: lum.map(l => Math.log10(l)),  // Logarithmic scale for Luminosity
                mode: 'markers',
                marker: {
                    // size: tmag.map(t => 20 - t),  // Increase the base size to make stars larger
                    // size: tmag.map(t => Math.max(5, 20 - t)),  // Minimum size of 5
                    size: tmag.map(t => 30 - t * 1.5),  // More aggressive scaling factor so that brighter stars appear significantly larger.
                    color: teff,  // Color based on temperature                
                    colorscale: [
                        [0, 'red'],
                        [0.5, 'yellow'],
                        [1, 'blue']
                    ],  // Blue to red (hot to cool stars)
                    opacity: 0.8,  // Higher opacity for better visibility
                    showscale: true,  // Enables the color scale bar
                    colorbar: {
                        title: {
                            text: 'Temperature (K)',
                            side: 'right'  // Position of the title
                        },
                        tickmode: 'array',  // Custom tick mode for specific temperature points
                        ticktext: ['3000 K', '5000 K', '7000 K', '10000 K', '30000 K']  // Labels for the temperatures ticks
                    },
                },
                // adding hover over text
                text: starData.map(star => `
                    <br><b>Star Basics</b>
                    <br>    Temperature: ${star.Teff} K
                    <br>    Star Type: ${star.c_star_type || 'Unknown'}
                    <br>    Luminosity Class: ${star.lumclass || 'Unknown'}
                    <br>    Luminosity: ${star.lum} L&#9737;
                    <br>    TESS magnitude: ${star.Tmag}
                    <br>    Exoplanet Host: ${star.c_has_exoplanet ? 'Yes' : 'No'}
                    <br><br><b>Comparison to the Sun</b>
                    <br>    Distance: ${star.distance || 'N/A'} pc
                    <br>    Mass: ${star.mass || 'N/A'} M&#9737;
                    <br>    Radius: ${star.radius || 'N/A'} R&#9737;
                    <br>    Metallicity: ${star.metallicity || 'N/A'} [Fe/H]
                    <br>
                `),
                hoverinfo: 'text'
            };

            const hrLayout = {
                title: 'Hertzsprung-Russell Diagram',
                xaxis: { title: 'Effective Temperature (K)', autorange: 'reversed' },  // Temperature decreasing from left to right
                yaxis: { title: 'Log Luminosity (L/Lsun)' },
                showlegend: false,
            };

            // Plot H-R diagram
            Plotly.newPlot('hr-plot', [hrTrace], hrLayout);

        })
        .catch(error => console.error('Error fetching data:', error));
});


/**
 * To display light curve data
 * Called by index.html, when user submits the "TIC ID" form: <button onclick="analyzeStar()">Analyze Star</button>
 * 
 * @returns - html block with data to display light curve info 
 */

function analyzeStar() {
    // Gets the TIC ID entered by the user and validate input not empty 
    const ticIdInput = 25155310 // document.getElementById('tic-id-input').value;
    if (!ticIdInput) {
        alert('Please enter a TIC ID.');
        return;
    }

    // Send POST request with TIC ID 
    fetch('/analyze_star', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ tic_id: ticIdInput })
    })
    .then(response => response.json())
    .then(data => {
        if (data.error) {
            alert(data.error);
            return;
        }

        // Access high-level info box data
        const metaObservationDetails = data.meta_observation_details;
        
        // Access high-level info box data
        const metaStellarDetails = data.meta_stellar_details;

        // Access detailed light curve data directly
        const lightCurve = data.light_curve;

        // Observation Details 
        document.getElementById('meta-observation-details').innerHTML = `
            <p><strong>Observation Start (BJD):</strong> ${metaObservationDetails.observation_start_bjd}</p>
            <p><strong>Observation Duration:</strong> ${metaObservationDetails.observation_duration} days (${metaObservationDetails.observation_start_date} - ${metaObservationDetails.observation_end_date})</p>
            <p><strong>Observation Cadence:</strong> ${metaObservationDetails.cadence_interval_minutes}-minute intervals</p>
        `;

        // Exoplanet Characteristics - Display the calculated properties with explanations
        document.getElementById('meta-exoplanet-details').innerHTML = `
            <p><strong>Orbital Period:</strong> ${lightCurve.period.toFixed(2)} days (Time taken to orbit around its star)</p>
            <p><strong>Transit Depth:</strong> ${(lightCurve.depth * 100).toFixed(2)}% (Decrease in star's brightness during transit)</p>
            <p><strong>Planet Radius:</strong> ${lightCurve.planet_radius.toFixed(2)} Solar Radii (Relative to Sun's size)</p>
        `;

        // Star Characteristics - Display properties with explanations
        document.getElementById('meta-stellar-details').innerHTML = `
            <p><strong>Radius:</strong> ${metaStellarDetails.rad.toFixed(2)} Solar Radii&#9737</p>
            <p><strong>Mass:</strong> ${metaStellarDetails.mass.toFixed(2)} Solar Masses&#9737</p>
            <p><strong>Distance:</strong> ${metaStellarDetails.d.toFixed(1)} parsecs</p>
            <p><strong>TESS Magnitude (Tmag):</strong> ${metaStellarDetails.Tmag.toFixed(2)}</p>
            <p><strong>Temperature (Teff):</strong> ${metaStellarDetails.Teff.toFixed(2)} K</p>
        `;
    
        /** 
         * ----- Start: Plotly setup - Updated light curve display using Plotly 
         *      Plotly Plot: Created a scatter plot with phase on the x-axis and normalized flux on the y-axis.
         *      Annotations: Added an annotation pointing to the transit dip.
         *      Interactivity: Users can zoom in and hover over data points to see exact values.
        */
        // Prepare data for Plotly
        const phase = lightCurve.phase;           // for folded light curve
        const flux_phase = lightCurve.flux_phase; // for folded light curve
        const time = lightCurve.time;             // for time-series light curve
        const flux_time = lightCurve.flux_time;  // for time-series light curve
        
        // Create the folded light curve plot
        const traceFolded = {
            x: phase,
            y: flux_phase,
            mode: 'markers',
            marker: {
                size: 5,
                color: 'blue'
            },
            name: 'Folded Light Curve'
        };

        const layoutFolded = {
            title: `Folded Light Curve for TIC ${ticIdInput} (RA:${metaStellarDetails.ra.toFixed(2)}, DEC:${metaStellarDetails.dec.toFixed(2)}, Sector:${lightCurve.sector})`,
            
            // phase on x axis
            xaxis: {
                title: 'Phase (Orbital Period Units)',
                showgrid: true,
                zeroline: true
            },
            
            // normalized flux on the y-axis
            yaxis: {
                title: 'Relative Brightness', // 'Normalized Flux',
                showgrid: true,
                zeroline: true
            },
            showlegend: false, // this turns on "light curve" toggle on the top right of the graph
            
            // annotation pointing to the transit dip
            annotations: [
                {
                    x: 0, // Adjust this based on where you typically observe dips
                    y: Math.min(...flux_phase),
                    xref: 'x',
                    yref: 'y',
                    text: 'Transit Dip',
                    showarrow: true,
                    arrowhead: 2,
                    ax: 90,
                    ay: 20
                }
            ],

            // width: 1300, // set based on how wide (wide U) we want the transit graph to be
            // height: 600, // set based on how long (deep U) we want the transit graph to be

        };

        const traceTimeSeries = {
            x: time,
            y: flux_time,
            mode: 'markers',
            marker: {
                size: 3,
                color: 'blue'
            },
            name: 'Time-Series Light Curve'
        };

        const layoutTimeSeries = {
            title: `Time-Series Light Curve for TIC ${ticIdInput} (RA:${metaStellarDetails.ra.toFixed(2)}, DEC:${metaStellarDetails.dec.toFixed(2)}, Sector:${lightCurve.sector})`,
            xaxis: {
                title: 'Time (Barycentric Julian Date - 2457000)',  // BJD stands for Barycentric Julian Date, a standard time system in astronomy.
                showgrid: true,
                zeroline: true
            },
            yaxis: {
                title: 'Normalized Flux',
                showgrid: true,
                zeroline: true
            },
            showlegend: false
        };

        // Plot both Folded and Time-Series Light Curves
        Plotly.newPlot('folded-light-curve-plot', [traceFolded], layoutFolded);
        Plotly.newPlot('light-curve-plot', [traceTimeSeries], layoutTimeSeries);

        // ----- End: Plotly setup

    })
    .catch(error => {
        console.error('Error:', error);
        alert('An error occurred while analyzing the star.');
    });
}