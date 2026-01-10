"""
Command-line interface for asimov-pastro.

This module provides CLI commands for computing p_astro and coherence
on gravitational wave posterior samples.
"""

import click
import json
import sys
import logging
from pathlib import Path


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--quiet', '-q', is_flag=True, help='Suppress all output except errors')
def cli(verbose, quiet):
    """
    asimov-pastro: Compute astrophysical probability for GW candidates.

    This tool computes p_astro (astrophysical probability) from posterior
    samples of gravitational wave parameter estimation.
    """
    # Setup logging
    if quiet:
        level = logging.ERROR
    elif verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )


@cli.command()
@click.argument('sample_file', type=click.Path(exists=True))
@click.option('--ifos', '-i', multiple=True,
              help='Interferometers to use (e.g., -i H1 -i L1). If not specified, will try to infer.')
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Output JSON file (default: <sample_file>_coherence.json)')
@click.option('--label', '-l', default=None,
              help='Label for PESummary files (e.g., C01:Mixed). If not specified, will try to auto-detect.')
def coherence(sample_file, ifos, output, label):
    """
    Compute signal coherence from posterior samples.

    This command calculates coherence metrics that distinguish astrophysical
    signals from instrumental noise by checking consistency across detectors.

    Example:

        asimov-pastro coherence posterior_samples.h5 -i H1 -i L1 -i V1
    """
    logger = logging.getLogger(__name__)

    # Convert to Path
    sample_path = Path(sample_file)

    # Generate output filename if not specified
    if output is None:
        output = sample_path.parent / f"{sample_path.stem}_coherence.json"
    else:
        output = Path(output)

    click.echo("=" * 70)
    click.echo("asimov-pastro: Coherence Calculation")
    click.echo("=" * 70)
    click.echo()
    click.echo(f"Input file:  {sample_path}")
    click.echo(f"Output file: {output}")
    click.echo()

    # Create calculator (import here to avoid circular imports)
    from .calculator import PastroCalculator
    calc = PastroCalculator(logger=logger)

    # Load samples
    click.echo("Step 1: Loading posterior samples")
    click.echo("-" * 70)

    try:
        # If label specified, we need to handle PESummary format specially
        if label:
            import h5py
            with h5py.File(sample_file, 'r') as f:
                if label not in f:
                    click.echo(f"Error: Label '{label}' not found in file", err=True)
                    click.echo(f"Available labels: {list(f.keys())}", err=True)
                    sys.exit(1)

                group = f[label]
                if 'posterior_samples' not in group:
                    click.echo(f"Error: No posterior_samples in {label}", err=True)
                    sys.exit(1)

                posterior = group['posterior_samples']
                samples = {}
                for field in posterior.dtype.names:
                    samples[field] = posterior[field][:]
        else:
            samples = calc.load_samples([str(sample_path)])

        click.echo(f"✓ Loaded {len(samples.get('ra', []))} posterior samples")
        click.echo(f"✓ Available parameters: {len(samples.keys())}")
        click.echo()

        # Show key parameters
        key_params = ['ra', 'dec', 'geocent_time', 'luminosity_distance',
                     'mass_1', 'mass_2', 'psi', 'theta_jn']
        missing_params = []
        for param in key_params:
            if param in samples:
                click.echo(f"  ✓ {param}")
            else:
                click.echo(f"  ✗ {param} (missing)")
                missing_params.append(param)

        if missing_params:
            click.echo()
            click.echo(f"Warning: Missing parameters: {', '.join(missing_params)}", err=True)
            click.echo("Coherence calculation may be limited.", err=True)

        click.echo()

    except Exception as e:
        click.echo(f"✗ Failed to load samples: {e}", err=True)
        if logger.level <= logging.DEBUG:
            import traceback
            traceback.print_exc()
        sys.exit(1)

    # Determine IFOs
    if ifos:
        ifos = list(ifos)
    else:
        # Try to infer from sample parameter names
        possible_ifos = ['H1', 'L1', 'V1', 'K1']
        ifos = []
        for ifo in possible_ifos:
            if any(ifo in key for key in samples.keys()):
                ifos.append(ifo)

        # Default to H1, L1 if can't infer
        if not ifos:
            ifos = ['H1', 'L1']
            click.echo(f"Warning: Could not infer IFOs, defaulting to {ifos}", err=True)

    # Compute coherence
    click.echo("Step 2: Computing coherence")
    click.echo("-" * 70)
    click.echo(f"Using interferometers: {', '.join(ifos)}")
    click.echo()

    try:
        coherence_results = calc.compute_coherence(samples, ifos=ifos)

        click.echo()
        click.echo("=" * 70)
        click.echo("RESULTS")
        click.echo("=" * 70)
        click.echo()

        # Display results with color
        def format_value(key, value):
            if value is None:
                return click.style("N/A", fg='yellow')
            elif isinstance(value, float):
                if 'consistency' in key or 'coherence' in key:
                    # Color code consistency values
                    if value > 0.8:
                        color = 'green'
                    elif value > 0.5:
                        color = 'yellow'
                    else:
                        color = 'red'
                    return click.style(f"{value:.4f}", fg=color, bold=True)
                else:
                    return f"{value:.4f}"
            else:
                return str(value)

        click.echo(f"Time delay consistency:  {format_value('consistency', coherence_results['time_delay_consistency'])}")
        click.echo(f"Amplitude consistency:   {format_value('consistency', coherence_results['amplitude_consistency'])}")

        if coherence_results['network_snr'] is not None:
            snr_color = 'green' if coherence_results['network_snr'] > 10 else 'yellow'
            click.echo(f"Network SNR (estimated): {click.style(f\"{coherence_results['network_snr']:.2f}\", fg=snr_color)}")
        else:
            click.echo("Network SNR:             N/A")

        click.echo(f"Overall coherence:       {format_value('coherence', coherence_results['overall_coherence'])}")
        click.echo()

        # Interpretation
        interp = coherence_results['interpretation']
        if 'HIGH' in interp:
            interp_color = 'green'
        elif 'MODERATE' in interp:
            interp_color = 'yellow'
        else:
            interp_color = 'red'

        click.echo(f"Interpretation: {click.style(interp, fg=interp_color, bold=True)}")
        click.echo()

        # Save results
        with open(output, 'w') as f:
            json.dump(coherence_results, f, indent=2)

        click.echo(f"✓ Results saved to: {click.style(str(output), fg='cyan')}")
        click.echo()

        # Exit with appropriate code
        if coherence_results['overall_coherence'] > 0.8:
            sys.exit(0)  # High coherence
        elif coherence_results['overall_coherence'] > 0.5:
            sys.exit(1)  # Moderate coherence
        else:
            sys.exit(2)  # Low coherence

    except Exception as e:
        click.echo(f"✗ Failed to compute coherence: {e}", err=True)
        if logger.level <= logging.DEBUG:
            import traceback
            traceback.print_exc()
        sys.exit(3)


@cli.command()
@click.argument('sample_file', type=click.Path(exists=True))
@click.option('--ifos', '-i', multiple=True,
              help='Interferometers to use (e.g., -i H1 -i L1)')
@click.option('--omicron', type=click.Path(exists=True), default=None,
              help='Omicron trigger file (optional)')
@click.option('--output', '-o', type=click.Path(), default=None,
              help='Output JSON file')
def pastro(sample_file, ifos, omicron, output):
    """
    Compute p_astro (astrophysical probability) from posterior samples.

    This combines coherence metrics with Bayesian odds to compute the
    probability that a candidate is astrophysical vs. terrestrial noise.

    Example:

        asimov-pastro pastro posterior_samples.h5 -i H1 -i L1
    """
    click.echo("=" * 70)
    click.echo("asimov-pastro: P_astro Calculation")
    click.echo("=" * 70)
    click.echo()
    click.echo(click.style("Note: Full p_astro calculation not yet implemented", fg='yellow'))
    click.echo("Currently only coherence is computed.")
    click.echo()
    click.echo("Run 'asimov-pastro coherence' for coherence calculation.")
    click.echo()
    sys.exit(1)


@cli.command()
@click.argument('sample_files', nargs=-1, type=click.Path(exists=True), required=True)
@click.option('--ifos', '-i', multiple=True,
              help='Interferometers to use for all files')
@click.option('--output-dir', '-o', type=click.Path(), default='.',
              help='Output directory for results (default: current directory)')
@click.option('--summary', type=click.Path(), default='coherence_summary.csv',
              help='Summary CSV file (default: coherence_summary.csv)')
def batch(sample_files, ifos, output_dir, summary):
    """
    Batch process multiple posterior sample files.

    Computes coherence for multiple events and generates a summary.

    Example:

        asimov-pastro batch GW*.h5 -i H1 -i L1 -o results/
    """
    import csv
    from datetime import datetime

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    click.echo("=" * 70)
    click.echo("asimov-pastro: Batch Processing")
    click.echo("=" * 70)
    click.echo()
    click.echo(f"Processing {len(sample_files)} files")
    click.echo(f"Output directory: {output_dir}")
    click.echo()

    results = []

    # Import here to avoid circular imports
    from .calculator import PastroCalculator

    # Process each file
    for i, sample_file in enumerate(sample_files, 1):
        sample_path = Path(sample_file)
        click.echo(f"[{i}/{len(sample_files)}] Processing {sample_path.name}...")

        output_file = output_dir / f"{sample_path.stem}_coherence.json"

        try:
            calc = PastroCalculator()
            samples = calc.load_samples([str(sample_path)])

            ifo_list = list(ifos) if ifos else ['H1', 'L1']
            coherence_result = calc.compute_coherence(samples, ifos=ifo_list)

            # Save individual result
            with open(output_file, 'w') as f:
                json.dump(coherence_result, f, indent=2)

            # Add to summary
            results.append({
                'file': sample_path.name,
                'n_samples': coherence_result['n_samples'],
                'ifos': ','.join(coherence_result['ifos']),
                'time_delay': coherence_result['time_delay_consistency'],
                'amplitude': coherence_result['amplitude_consistency'],
                'network_snr': coherence_result.get('network_snr', 'N/A'),
                'overall': coherence_result['overall_coherence'],
                'interpretation': coherence_result['interpretation'],
            })

            # Quick status
            if coherence_result['overall_coherence'] > 0.8:
                status = click.style('HIGH', fg='green')
            elif coherence_result['overall_coherence'] > 0.5:
                status = click.style('MODERATE', fg='yellow')
            else:
                status = click.style('LOW', fg='red')

            click.echo(f"  → Coherence: {status} ({coherence_result['overall_coherence']:.4f})")

        except Exception as e:
            click.echo(f"  ✗ Error: {e}", err=True)
            results.append({
                'file': sample_path.name,
                'n_samples': 'N/A',
                'ifos': 'N/A',
                'time_delay': 'N/A',
                'amplitude': 'N/A',
                'network_snr': 'N/A',
                'overall': 'N/A',
                'interpretation': f'ERROR: {e}',
            })

    click.echo()
    click.echo("=" * 70)
    click.echo("SUMMARY")
    click.echo("=" * 70)
    click.echo()

    # Write summary CSV
    summary_path = output_dir / summary
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'file', 'n_samples', 'ifos', 'time_delay', 'amplitude',
            'network_snr', 'overall', 'interpretation'
        ])
        writer.writeheader()
        writer.writerows(results)

    click.echo(f"Processed {len(sample_files)} files")
    click.echo(f"Summary saved to: {click.style(str(summary_path), fg='cyan')}")
    click.echo()

    # Quick statistics
    successful = [r for r in results if r['overall'] != 'N/A']
    if successful:
        high = sum(1 for r in successful if r['overall'] > 0.8)
        moderate = sum(1 for r in successful if 0.5 < r['overall'] <= 0.8)
        low = sum(1 for r in successful if r['overall'] <= 0.5)

        click.echo("Classification:")
        click.echo(f"  HIGH:     {high:3d} ({high/len(successful)*100:.1f}%)")
        click.echo(f"  MODERATE: {moderate:3d} ({moderate/len(successful)*100:.1f}%)")
        click.echo(f"  LOW:      {low:3d} ({low/len(successful)*100:.1f}%)")


@cli.command()
def version():
    """Show version information."""
    import asimov_pastro
    click.echo(f"asimov-pastro version: 0.1.0")
    click.echo(f"Python: {sys.version.split()[0]}")
    click.echo()
    click.echo("For more information: https://github.com/...")


def coherence_cmd():
    """Standalone coherence command to avoid circular imports."""
    # Import Click decorators inline
    import click

    @click.command()
    @click.argument('sample_file', type=click.Path(exists=True))
    @click.option('--ifos', '-i', multiple=True)
    @click.option('--output', '-o', type=click.Path(), default=None)
    @click.option('--label', '-l', default=None)
    @click.option('--verbose', '-v', is_flag=True)
    def run_coherence(sample_file, ifos, output, label, verbose):
        """Compute signal coherence from posterior samples."""
        # Just call the original function
        from click.testing import CliRunner
        ctx = click.Context(coherence)
        ctx.params = {
            'sample_file': sample_file,
            'ifos': ifos,
            'output': output,
            'label': label,
        }
        coherence.invoke(ctx)

    run_coherence()


if __name__ == '__main__':
    cli()
