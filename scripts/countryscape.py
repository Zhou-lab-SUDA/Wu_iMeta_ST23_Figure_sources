
#!/usr/bin/env python3
"""
countryscape.py - Calculate distances between countries using their names or codes.
"""

from geopy.geocoders import Nominatim
from geopy.distance import geodesic
import pycountry
import click

def get_country_code(country_input):
    """
    Convert country name or code to ISO 3166-1 alpha-2 country code
    """
    try:
        # If input is already a 2-letter code
        if len(country_input) == 2 and country_input.isalpha():
            country = pycountry.countries.get(alpha_2=country_input.upper())
        # If input is a 3-letter code
        elif len(country_input) == 3 and country_input.isalpha():
            country = pycountry.countries.get(alpha_3=country_input.upper())
        # If input is a country name
        else:
            country = pycountry.countries.get(name=country_input)
            if not country:
                # Try searching with fuzzy matching
                countries = pycountry.countries.search_fuzzy(country_input)
                country = countries[0] if countries else None

        if country:
            return country.alpha_2
        else:
            raise ValueError(f"Could not find country: {country_input}")
    except Exception as e:
        raise ValueError(f"Invalid country input: {country_input}. Error: {str(e)}")

def get_country_coordinates(country_input):
    """
    Get coordinates for a country using either name or code
    """
    try:
        country_code = get_country_code(country_input)
        geolocator = Nominatim(user_agent="countryscape")
        
        # Search using the standardized country code
        location = geolocator.geocode(country_code)
        
        if location:
            return (location.latitude, location.longitude)
        else:
            raise ValueError(f"Could not find coordinates for {country_input}")
    except Exception as e:
        raise ValueError(str(e))

def calculate_distance(country1, country2, unit='km'):
    """
    Calculate the distance between two countries
    """
    try:
        coords_1 = get_country_coordinates(country1)
        coords_2 = get_country_coordinates(country2)
        
        # Get full country names for display
        country1_name = pycountry.countries.get(alpha_2=get_country_code(country1)).name
        country2_name = pycountry.countries.get(alpha_2=get_country_code(country2)).name
        
        distance = geodesic(coords_1, coords_2).kilometers
        
        # Convert to miles if requested
        if unit.lower() in ['mi', 'miles']:
            distance = distance * 0.621371
            unit_name = "miles"
        else:
            unit_name = "kilometers"
            
        return distance, country1_name, country2_name, unit_name
    except ValueError as e:
        raise ValueError(str(e))

@click.command()
@click.argument('country1')
@click.argument('country2')
@click.option('--unit', '-u', type=click.Choice(['km', 'mi'], case_sensitive=False),
              default='km', help='Distance unit (km for kilometers, mi for miles)')
@click.option('--verbose', '-v', is_flag=True, help='Show additional information including coordinates')
def main(country1, country2, unit, verbose):
    """
    Countryscape - Calculate the distance between two countries.
    
    Countries can be specified using their full names, 2-letter codes (ISO 3166-1 alpha-2),
    or 3-letter codes (ISO 3166-1 alpha-3).

    Examples:
        countryscape "United States" Canada
        countryscape US CA --unit mi
        countryscape USA CAN --verbose
        countryscape "New Zealand" Australia --unit km --verbose
    """
    try:
        distance, country1_name, country2_name, unit_name = calculate_distance(country1, country2, unit)
        
        if verbose:
            coords1 = get_country_coordinates(country1)
            coords2 = get_country_coordinates(country2)
            click.echo(f"\nCountry 1: {country1_name}")
            click.echo(f"Coordinates: {coords1[0]:.4f}, {coords1[1]:.4f}")
            click.echo(f"\nCountry 2: {country2_name}")
            click.echo(f"Coordinates: {coords2[0]:.4f}, {coords2[1]:.4f}")
            click.echo("\nDistance:")
        
        click.echo(f"{country1_name} → {country2_name}: {distance:.2f} {unit_name}")
        
    except ValueError as e:
        click.echo(f"Error: {str(e)}", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"An unexpected error occurred: {str(e)}", err=True)
        raise click.Abort()

if __name__ == "__main__":
    main()

